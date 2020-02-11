#!/bin/bash

python statsChannel.py $1

gnuplot plot_umean.gnu
gnuplot plot_uvwrms.gnu

