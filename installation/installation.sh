#! /usr/bin/env bash
#!/bin/sh -l


#
# Description: Script to prepare environment to be able to run quaisar pipeline
#
# Usage: ./quaisar_installation.sh location_for
#
# Output location: No output created
#
# Modules required: None
#
# v1.0 (3/24/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

# Check that the 2 required softwares are installed (Python3 with biopython and Singularity)
command -v python3 >/dev/null 2>&1 && echo "Python 3 is installed"
command -v singularity
