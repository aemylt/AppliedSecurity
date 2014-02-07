#!/bin/bash

make debug
sh test.sh stage1
sh test.sh stage2
sh test.sh stage3
sh test.sh stage4
make clean
