#!/bin/bash

SRC=$(dirname "$(readlink -f "$0")")

SRC=$(readlink -f "$SRC/../TransTermHP/transterm_hp_v2.06")

$SRC/transterm -p $SRC/expterm.dat "$@"
