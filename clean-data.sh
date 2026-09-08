#!/bin/bash

Rscript clean-data.R -y 2017 -w 15
Rscript clean-data.R -y 2025 -w 15
Rscript clean-data.R -y 2017 -w 30
Rscript clean-data.R -y 2025 -w 30
Rscript clean-data.R -y 2017 -w 60
Rscript clean-data.R -y 2025 -w 60
