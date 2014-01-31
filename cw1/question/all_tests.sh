#!/bin/bash

make all
sh test.sh stage1
sh test.sh stage2
sh test.sh stage4
make clean
