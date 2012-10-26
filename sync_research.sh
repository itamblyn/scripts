#!/bin/bash
rsync -av --delete -e ssh ~/desktop/research/* itamblyn@192.168.1.5:~/research/
rsync -av --delete -e ssh ~/desktop/research/* itamblyn@chaffey.phys.dal.ca:~/research/
