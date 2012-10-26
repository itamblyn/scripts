#!/bin/bash

cd /Volumes/Backup/Backups.backupdb

rsync -av --delete telemachos itamblyn@chaffey.phys.dal.ca:backups/TELEMACHOS
