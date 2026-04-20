'''
    Particle definitions based on AliceO2 PhysicsConstants.h
    https://github.com/AliceO2Group/AliceO2/blob/dev/Common/Constants/include/CommonConstants/PhysicsConstants.h
'''

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class Particle:
    name: str
    mass: float          # GeV/c^2
    pdg: int
    label: str           # ROOT/LaTeX label
    pid_index: Optional[int] = None  # O2 PID index, if applicable


PARTICLES = {
    # --- Leptons ---
    'El':       Particle('El',   0.000510999,  11,        'e',                  pid_index=0),
    'Mu':       Particle('Mu',   0.1056584,    13,        '#mu',                pid_index=1),
    'NuE':      Particle('NuE',  0.0,          12,        '#nu_{e}'),
    'NuMu':     Particle('NuMu', 0.0,          14,        '#nu_{#mu}'),
    'NuTau':    Particle('NuTau',0.0,          16,        '#nu_{#tau}'),
    'Tau':      Particle('Tau',  1.77686,      15,        '#tau'),
    'Photon':   Particle('Photon',   0.0,          22,        '#gamma',             pid_index=10),

    # --- Light mesons ---
    'Pi0':       Particle('Pi0',   0.1349768,    111,       '#pi^{0}',            pid_index=9),
    'Pi':       Particle('Pi',   0.1395704,    211,       '#pi',                pid_index=2),
    'Eta':      Particle('Eta',  0.547862,     221,       '#eta'),
    'Omega782': Particle('Omega782', 0.78266,  223,       '#omega'),
    'EtaPrime': Particle('EtaPrime', 0.95778,  331,       "#eta'"),

    # --- Kaons ---
    'Ka':       Particle('Ka',   0.493677,     321,       'K',                  pid_index=3),
    'K0':       Particle('K0',   0.497611,     311,       'K^{0}',              pid_index=11),
    'K0Short':  Particle('K0Short', 0.497611,  310,       'K^{0}_{S}'),
    'K0Long':   Particle('K0Long',  0.497611,  130,       'K^{0}_{L}'),
    'K0Star892':    Particle('K0Star892',   0.89555,  313,  'K^{*0}(892)'),
    'KPlusStar892': Particle('KPlusStar892', 0.89167, 323, 'K^{*+}(892)'),
    'Phi':      Particle('Phi',  1.019461,     333,       '#phi'),

    # --- Light baryons ---
    'Pr':       Particle('Pr',   0.9382721,    2212,      'p',                  pid_index=4),
    'Ne':       Particle('Ne',   0.9395654,    2112,      'n'),

    # --- Strange baryons ---
    'Lambda':       Particle('Lambda',   1.115683,     3122,      '#Lambda',            pid_index=12),
    'Lambda1520': Particle('Lambda1520', 1.519, 3124,     '#Lambda(1520)'),
    'SigmaPlus':  Particle('SigmaPlus',  1.18937,  3222, '#Sigma^{+}'),
    'Sigma0':     Particle('Sigma0',     1.192642, 3212, '#Sigma^{0}'),
    'SigmaMinus': Particle('SigmaMinus', 1.197449, 3112, '#Sigma^{-}'),
    'Xi':       Particle('Xi',   1.32171,      3312,      '#Xi^{-}',            pid_index=15),
    'Xi0':      Particle('Xi0',  1.31486,      3322,      '#Xi^{0}'),
    'Om':       Particle('Om',   1.67245,      3334,      '#Omega^{-}',         pid_index=16),

    # --- Nuclei & hypernuclei ---
    'De':       Particle('De',   1.87561294257, 1000010020, 'd',                pid_index=5),
    'Tr':       Particle('Tr',   2.80892113298, 1000010030, '^{3}H',            pid_index=6),
    'He':       Particle('He',   2.80839160743, 1000020030, '^{3}He',           pid_index=7),
    'Al':       Particle('Al',   3.7273794066,  1000020040, '^{4}He',           pid_index=8),
    'Li4':      Particle('Li4',  3.7513,        1000030040, '^{4}Li'),
    'H3L':       Particle('H3L',   2.991134,      1010010030, '^{3}_{#Lambda}H',  pid_index=13),
    'H4L':      Particle('H4L',  3.922434,      1010010040, '^{4}_{#Lambda}H',  pid_index=14),
    'HyperHe4': Particle('HyperHe4', 3.921728,  1010020040, '^{4}_{#Lambda}He'),
    'HyperHe5': Particle('HyperHe5', 4.839961,  1010020050, '^{5}_{#Lambda}He'),
    'HyperHe4Sigma': Particle('HyperHe4Sigma', 3.995, 1110020040, '^{4}_{#Sigma}He'),

    # --- Charm mesons ---
    'D0':       Particle('D0',    1.86484, 421,   'D^{0}'),
    'DPlus':    Particle('DPlus', 1.86966, 411,   'D^{+}'),
    'DS':       Particle('DS',    1.96835, 431,   'D^{+}_{s}'),
    'DSStar':   Particle('DSStar',2.1122,  433,   'D^{*+}_{s}'),
    'DStar':    Particle('DStar', 2.01026, 413,   'D^{*+}'),
    'DStar0':   Particle('DStar0',2.00685, 423,   'D^{*0}'),

    # --- Charm baryons ---
    'LambdaCPlus': Particle('LambdaCPlus', 2.28646, 4122, '#Lambda^{+}_{c}'),
    'XiCPlus':     Particle('XiCPlus',    2.46771, 4232, '#Xi^{+}_{c}'),
    'XiC0':        Particle('XiC0',       2.47044, 4132, '#Xi^{0}_{c}'),
    'OmegaC0':     Particle('OmegaC0',    2.6952,  4332, '#Omega^{0}_{c}'),
    'SigmaC0':     Particle('SigmaC0',    2.45375, 4112, '#Sigma^{0}_{c}'),
    'SigmaCPlusPlus': Particle('SigmaCPlusPlus', 2.45397, 4222, '#Sigma^{++}_{c}'),
    'XiCCPlusPlus':   Particle('XiCCPlusPlus',   3.62155, 4422, '#Xi^{++}_{cc}'),

    # --- Charmonia ---
    'JPsi':     Particle('JPsi',   3.0969,   443,      'J/#psi'),
    'ChiC1':    Particle('ChiC1',  3.51067,  20443,    '#chi_{c1}'),
    'X3872':    Particle('X3872',  3.87165,  9920443,  'X(3872)'),

    # --- Beauty mesons ---
    'BPlus':    Particle('BPlus',  5.27934,  521,  'B^{+}'),
    'B0':       Particle('B0',     5.27966,  511,  'B^{0}'),
    'BS':       Particle('BS',     5.36692,  531,  'B^{0}_{s}'),
    'BCPlus':   Particle('BCPlus', 6.27447,  541,  'B^{+}_{c}'),

    # --- Beauty baryons ---
    'LambdaB0': Particle('LambdaB0', 5.6196, 5122, '#Lambda^{0}_{b}'),
    'XiB0':     Particle('XiB0',     5.7919, 5232, '#Xi^{0}_{b}'),
}