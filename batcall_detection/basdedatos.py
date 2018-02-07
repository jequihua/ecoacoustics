# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 21:31:29 2018

@author: ttonaru
"""

import json
with open(".\muestra\metadatos_tona.json", "r") as f:
    datos = json.load(f)
    
print(datos[0]["_id"])

