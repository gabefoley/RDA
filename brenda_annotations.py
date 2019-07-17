annot_dict = {"APPLICATION": "AP", "ACTIVATING_COMPOUND": "AC", "CLONED": "CL", "COFACTOR": "CF", "ENGINEERING": "EN",
              "EXPRESSION" : "EXP", "GENERAL_STABILITY": "GS", "IC50_VALUE" : "IC50", "INHIBITORS": "IN",
              "KI_VALUE": "KI", "KM_VALUE": "KM" , "LOCALIZATION": "LO", "METALS_IONS": "ME", "MOLECULAR_WEIGHT": "MW",
              "NATURAL_SUBSTRATE_PRODUCT" : "NSP", "OXIDATION_STABILITY": "OS", "PH_OPTIMUM" : "PHO", "PH_RANGE" :
                  "PHR", "PURIFICATION": "PU", "REACTION": "RE", "REACTION_TYPE": "RT", "RECOMMENDED_NAME": "RN",
              "REFERENCE": "RF", "SOURCE_TISSUE": "ST", "SPECIFIC_ACTIVITY": "SA", "STORAGE_STABILITY": "SS",
              "SUBSTRATE_PRODUCT": "SP", "SUBUNITS": "SU", "SYNONYMS": "SY", "SYSTEMATIC_NAME": "SN",
              "TEMPERATURE_OPTIMUM": "TO", "TEMPERATURE_RANGE": "TR", "TEMPERATURE_STABILITY": "TS", "TURNOVER_NUMBER":
                  "TN"}


def get_brenda_dict(annot_list):
    """
    Create a cut down version of the annotation dictionary, based on what the user has actually asked for
    :param annot_list: List of annotations user is interested in
    :return: Dictionary mapping annotations to abbreviated name
    """
    brenda_dict = {}

    for annot in annot_list:
        brenda_dict[annot] = annot_dict[annot]

    return brenda_dict