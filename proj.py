import json
import math

## GLOBAL OPTS
DEBUG = False

## OBJECTS ##
class obj:
    def __init__(self, dct):
        self.__dict__.update(dct)

class section:
    def __init__(self, G, F, L, ASE, desc):
        self.G = G
        self.F = F
        self.L = L
        self.ASE = ASE
        self.desc = desc 

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

## FUNCTIONS ##
def conv_to_objects (config):
    cfg_obj = json.loads(json.dumps(config), object_hook=obj)
    ba = cfg_obj.BA
    a = cfg_obj.A
    pa = cfg_obj.PA
    ra = cfg_obj.RA
    fra = cfg_obj.FRA
    link = cfg_obj.link
    losses = cfg_obj.losses
    otu = cfg_obj.otu
    return ba, a, pa, ra, fra, link, losses, otu

def loss(losses, l):
    return ( losses.loss_per_km + losses.loss_margin_per_km) * l + losses.loss_margin_db 

def dispersion(link):
    l1_comp = l2_comp = lt_comp = 0
    while not (link.min_disp <= (link.l1 - l1_comp) * link.disp_per_km <= link.max_disp):
        l1_comp += 20
    while not (link.min_disp <= (link.l2 - l2_comp) * link.disp_per_km <= link.max_disp):
        l2_comp += 20
    while not (link.min_disp <= (link.l1 + link.l2 - lt_comp) * link.disp_per_km <= link.max_disp):
        lt_comp += 20
    return l1_comp, l2_comp, lt_comp 


def calculate_losses(p_max_a, p_max_ba, p_min_no_ra, p_min_ra, fra_gain):
    return [
        p_max_a - p_min_no_ra,
        p_max_ba - p_min_no_ra,
        p_max_a - p_min_ra,
        p_max_ba - p_min_ra,
        p_max_a - p_min_ra + fra_gain
    ]


def get_amplifiers(L1, L2, max_losses):
    i1, i2 = 0, 0
    configs = [
        "A - /",
        "BA - /",
        "A - RA",
        "BA - RA",
        "A/FRA - RA"

    ]
    for loss in max_losses:
        if L1 <= loss:
            break
        i1 += 1
    for loss in max_losses:
        if L2 <= loss:
            break
        i2 += 1

    print("Line 1: " + bcolors.OKCYAN + configs[i1] + bcolors.ENDC + "\nLine 2: " + bcolors.OKCYAN + configs[i2] + bcolors.ENDC + "\n") 
    return [
            (i1 % 2, 0 if i1 < 2 else 1, 1 if i1 == 4 else 0),
            (i2 % 2, 0 if i2 < 2 else 1, 1 if i2 == 4 else 0)
        ] # [A/AB, -/RA, -/FRA]

def get_voa_mcd_attenuation(config, L, p_max_a, p_max_ba, a, ba, pa, ra, fra, losses, p_boost_out):
    p_tx = p_max_a if config[0] == 0 else p_max_ba
    p_rx = p_tx - L + (ra.gain if config[1] else 0) + (fra.gain if config[2] else 0)
    p_voa = p_rx - losses.lp_rx - losses.fiu + pa.gain
    p_boost_in = p_boost_out - a.gain 
    voa_mcd_att = p_voa - p_boost_in 
    return voa_mcd_att

def get_dispersion_modules(att1, att2, link, losses):
    ## we are absolutely brute forcing this XD
    flag_precomp_1 = flag_precomp_2 = False
    mcd1_opts = []
    mcd2_opts = []
    mcd_losses = [0, losses.comp_20, losses.comp_40, losses.comp_60, losses.comp_80]
    for i in range(2): 
        # calculate possible modules for MCD\n2
        mcd2 = 0
        while not (link.min_disp <= (link.l2 - mcd2 - 20 * i) * link.disp_per_km <= link.max_disp):
            mcd2 += 20
        while (link.min_disp <= (link.l2 - mcd2 - 20 * i ) * link.disp_per_km <= link.max_disp):
            mcd_att = mcd_losses[mcd2 // 20 % 5]
            mcd_att += losses.comp_100 * (mcd2 // 100)
            if att2 > mcd_att:  
                mcd2_opts.append((mcd2, (link.l2 - mcd2 - 20 * i) * link.disp_per_km, i == 1)) # last flag denotes precompensation
            mcd2 += 20

        # calculate possible modules for MCD1
        mcd1 = 0
        while not (link.min_disp <= (link.l1 - mcd1 - 20 * i) * link.disp_per_km <= link.max_disp):
            mcd1 += 20
        while (link.min_disp <= (link.l1 - mcd1 - 20 * i) * link.disp_per_km <= link.max_disp):
            mcd_att = mcd_losses[mcd1 // 20 % 5]
            mcd_att += losses.comp_100 * (mcd1 // 100)
            if att1 > mcd_att:  
                mcd1_opts.append((mcd1, (link.l1 - mcd1 - 20 * i) * link.disp_per_km, i == 1))
            mcd1 += 20

    # calculate the possible combinations MCD1, MCD2 that will work
    mcd_opts = []
    for mcd1, disp1, flag_precomp_1 in mcd1_opts:
        for mcd2, disp2, flag_precomp_2 in mcd2_opts:
            if link.min_disp < disp1 + disp2 < link.max_disp:
                    mcd_opts.append((mcd1, mcd2, disp1, disp2, flag_precomp_1, flag_precomp_2))
   
    return mcd_opts

def get_osnr(p_ch_a, p_ch_ba, L1, L2, amps_l1, amps_l2, losses, voa_mcd_1, voa_mcd_2, a, ba, pa, ra, fra, flag_pc_OADM):
    NO_ASE = -999999 # very low dBm value for ASE noise
    ASE_LIM = -1000 # limit to print NO_ASE

    # amps == [A/AB, -/RA, -/FRA]
    sections = [] # G, F, L, ASE
    print(f"{bcolors.BOLD}OTM1: PA - {"A" if amps_l1[0] == 0 else "BA"} - {"/" if amps_l1[2] == 0 else "FRA"}{bcolors.ENDC}")
    # OTM 1
    # OA-MCD
    sections.append(section(pa.gain, pa.f, pa.gain, pa.ase, "PA + VOA")) # PA 
    booster, b_desc = (a, "A + IU + LP") if amps_l1[0] == 0 else (ba, "BA + IU + LP")
    sections.append(section(booster.gain, booster.f, losses.fiu + losses.lp_tx, booster.ase, b_desc)) # BOOSTER
    # FRA
    sections.append(section(0, 0, L2, NO_ASE, "Fiber 1") if amps_l1[2] == 0 else section(fra.gain, fra.f, L1, fra.ase, "FRA + Fiber 1"))

    print(f"{bcolors.BOLD}OADM({"OBA" if not flag_pc_OADM else "MCA-CD"}): {"/" if amps_l1[1] == 0 else "RA"} - PA - A - MUX|DEMUX - {"/" if not flag_pc_OADM else "PA"} - {"A" if amps_l2[0] == 0 else "BA"} - {"/" if amps_l2[2] == 0 else "FRA"}{bcolors.ENDC}")
    # OADM 
    # RA
    sections.append(section(0, 0, losses.lp_rx + losses.fiu, NO_ASE,"LP + IU") if amps_l1[1] == 0 else section(ra.gain, ra.f, losses.lp_rx + losses.fiu, ra.ase, "RA + LP + IU"))
    # OA-MCD
    sections.append(section(pa.gain, pa.f, voa_mcd_1, pa.ase, "PA + VOA")) 
    booster, b_desc = a, "A + DEMUX + MUX + VOA"
    sections.append(section(booster.gain, booster.f, p_ch_a - pa.min_in, booster.ase, b_desc))
    
    booster, b_desc = (a, "A") if amps_l2[0] == 0 else (ba, "BA")
    if flag_pc_OADM:
        # PA if using OA-MCD
        sections.append(section(pa.gain, pa.f, pa.gain, pa.ase, "PA + VOA")) 
    sections.append(section(booster.gain, booster.f, losses.fiu + losses.lp_tx, booster.ase, b_desc))
    # FRA 
    sections.append(section(0, 0, L2, NO_ASE, "Fiber 2") if amps_l2[2] == 0 else section(fra.gain, fra.f, L2, fra.ase, "FRA + Fiber 2"))

    print(f"{bcolors.BOLD}OTM2: {"/" if amps_l2[1] == 0 else "RA"} - A{bcolors.ENDC}\n")
    # OTM 2
    # RA
    sections.append(section(0, 0, losses.lp_rx + losses.fiu, NO_ASE, "LP + IU") if amps_l2[1] == 0 else section(ra.gain, ra.f, losses.lp_rx + losses.fiu, ra.ase, "RA + LP + IU"))
    # OA-MCD
    sections.append(section(pa.gain, pa.f, voa_mcd_2, pa.ase, "PA + VOA")) # PA 
    booster = a
    booster, b_desc = a, "A"
    sections.append(section(booster.gain, booster.f, 0, booster.ase, b_desc)) # BOOSTER, last section has L = 0

    ase = []
    if DEBUG:
        print("i, G, L, ASE")
        for i, s in enumerate(sections):
            print(f"{bcolors.WARNING}{i}: {s.G:.3f}, {s.L:.3f}, {s.ASE:.3f}, {bcolors.ENDC}")

    for sect in reversed(sections):
        ase.append((ase[-1] if ase else 0) + sect.G - sect.L)
    ase.reverse()

    print(f"{bcolors.UNDERLINE}P_ase per section is:{bcolors.ENDC}")
    for i, sect in enumerate(sections):
        # gain was added to make the algo faster, now is being subtracted
        # we are also converting ASE back from logarithm
        print(f"{bcolors.OKCYAN}P_ase_{i} = {ase[i] - sect.G + sect.ASE:.2f}dBm{bcolors.ENDC}, {sect.desc}" if sect.ASE > ASE_LIM 
              else f"{bcolors.OKCYAN}P_ase_{i} = NONE{bcolors.ENDC}")

        ase[i] = 10 ** (0.1 * (ase[i] - sect.G + sect.ASE)) 

    ase_tot = 10 * math.log10(sum(ase)) # total ASE in dB

    osnr = p_ch_a - ase_tot
    print(f"Total ASE of the system is:\n{bcolors.OKBLUE}P_ase = {ase_tot:.2f}dB{bcolors.ENDC}")
    print(f"OSNR of the system is:\n{bcolors.OKBLUE}OSNR = {osnr:.2f}dB{bcolors.ENDC}")

def print_opts(mcd_opts, print_all = True):
    if not mcd_opts:
        print(bcolors.FAIL + "NO OPTIONS FOUND\n")
        return
    print(f"Following options are available: ({bcolors.OKGREEN} valid, {bcolors.WARNING} if no green, {bcolors.FAIL} if no yellow {bcolors.ENDC})")
    flag_printed = False
    for c in [bcolors.OKGREEN, bcolors.WARNING, bcolors.FAIL]:
        if flag_printed and not print_all:
            break
        for opt in mcd_opts:
            mcd1, mcd2, disp1, disp2, flag_precomp_1, flag_precomp_2 = opt 
            color = bcolors.FAIL if flag_precomp_1 and flag_precomp_2 else bcolors.WARNING if flag_precomp_1 or flag_precomp_2 else bcolors.OKGREEN
            # this ensures we print in the order GREEN -> YELLOW -> RED
            if c != color:
                continue
            mcd1_1 = mcd1 if mcd1 <= 100 else 100
            mcd1_2 = 0 if mcd1 <= 100 else mcd1- 100 
            mcd2_1 = mcd2 if mcd2 <= 100 else 100 
            mcd2_2 = 0 if mcd2 <= 100 else mcd2 - 100 
            print(f"""{color}MCD{mcd1_1}{f"+{mcd1_2}" if mcd1_2 else ""} and MCD{mcd2_1}{f"+{mcd2_2}" if mcd2_2 else ""}
            link1 dispersion = {disp1}ps/nm
            link2 dispersion = {disp2}ps/nm
            total dispersion = {disp1 + disp2}ps/nm
            precompensation: {"None" if not flag_precomp_1 and not flag_precomp_2 else ""}{"link 1" if flag_precomp_1 else ""}{" and " if flag_precomp_1 and flag_precomp_2 else ""}{"link 2" if flag_precomp_2 else ""}{bcolors.ENDC}\n""")
            flag_printed = True


## TEST FUNCTIONS ##
def test_loss(losses):
    # 2
    losses.loss_per_km = 0.271 
    assert(round(loss(losses, 59),2) == 20.17)
    # 5
    losses.loss_per_km = 0.231
    assert(round(loss(losses, 108),2) == 30.11)
    # 14
    losses.loss_per_km = 0.243
    assert(round(loss(losses, 148), 2) == 41.92)

    print(f"{bcolors.OKGREEN}LOSS TESTS PASSED{bcolors.ENDC}")

def main(TEST):
    # load the project configuration data 
    config = {}

    with open('conf.json', 'r') as file:
        config = json.load(file)

    ba, a, pa, ra, fra, link, losses, otu = conv_to_objects(config)

    ## call tests if in test mode #TODO : add more tests
    if TEST:
        test_loss(losses)
        return
    
    # calculate losses
    L1 = loss(losses, link.l1)
    L2 = loss(losses, link.l2)

    print("Line losses:\n{}L1 = {:.2f}, L2 = {:.2f}\n{}".format(bcolors.OKCYAN, L1, L2, bcolors.ENDC))

    # find max power per channel
    p_ch_ba = round( 10 * math.log10( 10**( ba.max_out / 10) / (otu.channel_count // 2)))
    p_ch_a = round( 10 * math.log10( 10**( a.max_out / 10) / (otu.channel_count // 2)))
    print("Max power per channel:\n{}BA = {}dB, A = {}dBm{}\n".format(bcolors.OKCYAN, p_ch_ba, p_ch_a, bcolors.ENDC))

    p_max_ba = p_ch_ba - losses.fiu - losses.lp_tx
    p_max_a = p_ch_a - losses.fiu - losses.lp_tx
    print("Max output power:\n{}BA = {}dBm, A = {}dBm{}\n".format(bcolors.OKCYAN, p_max_ba, p_max_a, bcolors.ENDC))

    p_min_no_ra = pa.min_in + losses.fiu + losses.lp_rx
    p_min_ra = p_min_no_ra - ra.gain
    print("Min reciever power:\n{}No RA = {}dBm, With RA = {}dBm{}\n".format(bcolors.OKCYAN, p_min_no_ra, p_min_ra, bcolors.ENDC))

    max_losses = calculate_losses(p_max_a, p_max_ba, p_min_no_ra, p_min_ra, fra.gain)
    print("Maximum line losses for each of the following combinations is:\n",
        f"""{bcolors.OKCYAN}
        A - / = {max_losses[0]}dB
        BA - / = {max_losses[1]}dB
        A - RA = {max_losses[2]}dB
        BA - RA = {max_losses[3]}dB
        A/FRA - RA = {max_losses[4]}dB
        {bcolors.ENDC}""")
    
    amps_l1, amps_l2 = get_amplifiers(L1, L2, max_losses)
    att1 = get_voa_mcd_attenuation(amps_l1, L1, p_max_a, p_max_ba, a, ba, pa, ra, fra, losses, p_ch_a)  
    att2 = get_voa_mcd_attenuation(amps_l2, L2, p_max_a, p_max_ba, a, ba, pa, ra, fra, losses, p_ch_a)  

    print(f"Attenuation on each of VOA + MCD is:{bcolors.OKCYAN}\nVOA_OADM={att1:.2f}\nVOA_OTM2={att2:.2f}\n{bcolors.ENDC}")
    mcd_opts = get_dispersion_modules(att1, att2, link, losses)
    print_opts(mcd_opts, print_all=False)

    # see if there is no options without pc2
    flag = all(pc2 for _,_,_,_,_,pc2 in mcd_opts) # only true if all options have precomp on 2
    get_osnr(p_ch_a, p_ch_ba, L1, L2, amps_l1, amps_l2, losses, att1, att2, a, ba, pa, ra, fra, flag_pc_OADM=flag)

if __name__ == "__main__":
    TEST = False
    main(TEST)
