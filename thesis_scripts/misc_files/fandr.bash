#!/bin/bash

# Assign the filename
filename=$1

# Take the search string
read -p "Enter the search string: " search

# Take the replace string
read -p "Enter the replace string: " replace

if [[ $search != "" && $replace != "" ]]; then
  sed -i '' "s/$search/$replace/" $filename
fi
