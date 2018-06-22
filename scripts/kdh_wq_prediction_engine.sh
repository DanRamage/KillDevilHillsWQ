#!/bin/bash

source /usr/local/virtualenv/pyenv-2.7.11/bin/activate

cd /home/xeniaprod/scripts/KillDevilHillsWQ/scripts;

python kdh_prediction_engine.py --ConfigFile=/home/xeniaprod/scripts/KillDevilHillsWQ/config/kdh_prediction_engine.ini >> /home/xeniaprod/tmp/log/kdh_wq_prediction_engine_sh.log 2>&1
