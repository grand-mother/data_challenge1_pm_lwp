import numpy as np
import electronic_chain.ec_config as ec_config
# from electronic_chain.XDU_electronic_chain.functions import *

def voltage2adc(traces_t, **kwargs):
    """Converts voltage traces to ADC traces"""

    Length_AD = 14  # Significant bit of the value, plus a sign bit in addition to this
    Range_AD = 1.8 * 1e6  # Vpp,unit:uv
    delta = Range_AD / 2 / (2 ** (Length_AD - 1))
    traces_t = np.sign(traces_t) * np.floor(abs(traces_t) / delta) * delta

    # Inside pipeline return - a dictionary
    if ec_config.in_pipeline:
        return {"traces_t": traces_t}
    # Outside pipeline return - raw values
    else:
        return traces_t

