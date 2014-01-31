#!/bin/bash

./modmul $1 < $1.input > output.txt
diff output.txt $1.output
rm output.txt
