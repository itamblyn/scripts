#!/bin/bash

wget http://automation.whatismyip.com/n09230945.asp

mv n09230945.asp telemachos

scp telemachos itamblyn@ndorrance.dyndns.org:machines
