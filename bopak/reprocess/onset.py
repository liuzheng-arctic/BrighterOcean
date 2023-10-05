import numpy as np
import pandas as pd

from .sic import sic

class onset(sic):
    
    def _time(self):
        return pd.to_datetime(self.fn.name[:4])