"""Prespecified study labels and comparison families."""
CENTRES = ("Main centre", "Centre I", "Centre II", "Centre III", "Centre IV")
EXTERNAL_N = {"Centre I": 14, "Centre II": 21, "Centre III": 14, "Centre IV": 25}
METRICS = ("DSC", "HD95", "ASD")
TARGETS = ("GTVp", "CTV1", "CTV2", "GTVn")
SMALL = ("Optic chiasm", "Cochlea", "Eustachian tube bone", "Internal auditory canal",
         "Lens", "Optic nerve", "Pituitary", "Temporomandibular joint", "Tympanic cavity",
         "Vestibule/semicircular canals")
LARGE = ("Eye", "Glottic larynx", "Supraglottic larynx", "Mandible", "Oral cavity",
         "Parotid gland", "Brain stem", "Thyroid", "Temporal lobe", "Submandibular gland", "Spinal cord")
GROUPS = {"Small OARs": SMALL, "Large OARs": LARGE}
GROUP_LABELS = {"Small OARs": "Small OARs (≤ 4 cm³)", "Large OARs": "Large OARs (> 4 cm³)"}
ENDPOINTS = (*GROUP_LABELS.values(), *TARGETS)
STRATEGIES = ("Image-only baseline", "Direct DIR propagation", "PCG-UNet")
TIMES = ("Algorithm-processing time (s)", "Review-and-modification interval (s)", "Derived combined time (s)")
CASE_ORDER = (
    ("Central nervous", "Brain stem"), ("Central nervous", "Spinal cord"),
    ("Central nervous", "Eye"), ("Central nervous", "Temporal lobe"),
    ("Visual", "Optic chiasm"), ("Visual", "Lens"), ("Visual", "Pituitary"), ("Visual", "Optic nerve"),
    ("Auditory", "Cochlea"), ("Auditory", "Eustachian tube bone"),
    ("Auditory", "Internal auditory canal"), ("Auditory", "Tympanic cavity"),
    ("Auditory", "Vestibule/semicircular canals"),
    ("Oral", "Parotid gland"), ("Oral", "Temporomandibular joint"), ("Oral", "Oral cavity"),
    ("Oral", "Submandibular gland"), ("Oral", "Mandible"),
    ("Pharyngo-laryngeal", "Supraglottic larynx"), ("Pharyngo-laryngeal", "Thyroid"),
    ("Pharyngo-laryngeal", "Glottic larynx"),
)


def models(centre):
    return (("Conventional U-Net comparator", "Locked source model") if centre == "Main centre"
            else ("Pre-adaptation model", "Locally adapted model"))


def expected_n(centre, roi):
    return (12 if roi in ("GTVp", "GTVn") else 29) if centre == "Main centre" else EXTERNAL_N[centre]

