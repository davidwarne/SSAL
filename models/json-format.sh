#!/bin/bash
FILE=$1
python -m json.tool $FILE > temp
mv temp $FILE

