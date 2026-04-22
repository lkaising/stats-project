"""Locked parameter values and design constants."""

from dataclasses import dataclass
from types import MappingProxyType
from typing import Literal, Mapping

SUBJECT_ID = "S01"
N_TRIALS = 12
SEED = 42

Site = Literal["dorsal_hand_L", "dorsal_hand_R", "antecubital_L", "antecubital_R"]
Wavelength = Literal[850, 940, 1050]
SingleBandCondition = Literal["850", "940", "1050"]
RatioCondition = Literal["850_940", "940_1050"]
Condition = SingleBandCondition | RatioCondition

SITES: tuple[Site, ...] = (
    "dorsal_hand_L",
    "dorsal_hand_R",
    "antecubital_L",
    "antecubital_R",
)
WAVELENGTHS: tuple[Wavelength, ...] = (850, 940, 1050)
SINGLE_BAND_CONDITIONS: tuple[SingleBandCondition, ...] = ("850", "940", "1050")
RATIO_CONDITIONS: tuple[RatioCondition, ...] = ("850_940", "940_1050")
CONDITIONS: tuple[Condition, ...] = SINGLE_BAND_CONDITIONS + RATIO_CONDITIONS

RATIO_PARENTS: dict[RatioCondition, tuple[SingleBandCondition, SingleBandCondition]] = {
    "850_940": ("850", "940"),
    "940_1050": ("940", "1050"),
}

MU_B: dict[SingleBandCondition, float] = {"850": 0.620, "940": 0.590, "1050": 0.550}
MU_V: dict[SingleBandCondition, float] = {"850": 0.508, "940": 0.519, "1050": 0.506}

SIGMA_V = 0.008
SIGMA_B = 0.007
SIGMA_R = 0.006

SITE_OFFSETS: dict[Site, float] = {
    "dorsal_hand_L": -0.012,
    "dorsal_hand_R": -0.010,
    "antecubital_L": +0.012,
    "antecubital_R": +0.010,
}

INTERACTIONS: dict[tuple[Site, Condition], float] = {
    ("dorsal_hand_L", "850"): -0.006,
    ("dorsal_hand_L", "940"): 0.000,
    ("dorsal_hand_L", "1050"): +0.006,
    ("dorsal_hand_L", "850_940"): -0.003,
    ("dorsal_hand_L", "940_1050"): +0.003,
    ("dorsal_hand_R", "850"): -0.005,
    ("dorsal_hand_R", "940"): +0.001,
    ("dorsal_hand_R", "1050"): +0.005,
    ("dorsal_hand_R", "850_940"): -0.002,
    ("dorsal_hand_R", "940_1050"): +0.004,
    ("antecubital_L", "850"): +0.006,
    ("antecubital_L", "940"): 0.000,
    ("antecubital_L", "1050"): -0.006,
    ("antecubital_L", "850_940"): +0.003,
    ("antecubital_L", "940_1050"): -0.003,
    ("antecubital_R", "850"): +0.005,
    ("antecubital_R", "940"): -0.001,
    ("antecubital_R", "1050"): -0.005,
    ("antecubital_R", "850_940"): +0.002,
    ("antecubital_R", "940_1050"): -0.004,
}


@dataclass(frozen=True)
class Parameters:
    subject_id: str
    n_trials: int
    seed: int
    sites: tuple[Site, ...]
    wavelengths: tuple[Wavelength, ...]
    single_band_conditions: tuple[SingleBandCondition, ...]
    ratio_conditions: tuple[RatioCondition, ...]
    conditions: tuple[Condition, ...]
    ratio_parents: Mapping[RatioCondition, tuple[SingleBandCondition, SingleBandCondition]]
    mu_B: Mapping[SingleBandCondition, float]
    mu_V: Mapping[SingleBandCondition, float]
    sigma_V: float
    sigma_B: float
    sigma_R: float
    site_offsets: Mapping[Site, float]
    interactions: Mapping[tuple[Site, Condition], float]


def load_parameters(seed: int = SEED) -> Parameters:
    return Parameters(
        subject_id=SUBJECT_ID,
        n_trials=N_TRIALS,
        seed=seed,
        sites=SITES,
        wavelengths=WAVELENGTHS,
        single_band_conditions=SINGLE_BAND_CONDITIONS,
        ratio_conditions=RATIO_CONDITIONS,
        conditions=CONDITIONS,
        ratio_parents=MappingProxyType(RATIO_PARENTS),
        mu_B=MappingProxyType(MU_B),
        mu_V=MappingProxyType(MU_V),
        sigma_V=SIGMA_V,
        sigma_B=SIGMA_B,
        sigma_R=SIGMA_R,
        site_offsets=MappingProxyType(SITE_OFFSETS),
        interactions=MappingProxyType(INTERACTIONS),
    )
