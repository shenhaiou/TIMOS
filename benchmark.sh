#!/bin/sh
time -p ./source/timos -f example/onelayer/one_layer_18_18_1_2.mesh -s example/onelayer/one_layer_18_18_1_2.source -p example/onelayer/mua005_mus05.opt -m is  -t 32
time -p ./source/timos -f example/onelayer/one_layer_18_18_1_2.mesh -s example/onelayer/one_layer_18_18_1_2.source -p example/onelayer/mua005_mus10.opt -m is  -t 32
time -p ./source/timos -f example/onelayer/one_layer_18_18_1_2.mesh -s example/onelayer/one_layer_18_18_1_2.source -p example/onelayer/mua020_mus05.opt -m is  -t 32
time -p ./source/timos -f example/onelayer/one_layer_18_18_1_2.mesh -s example/onelayer/one_layer_18_18_1_2.source -p example/onelayer/mua020_mus10.opt -m is  -t 32

time -p ./source/timos -f example/fourlayer/FourLayer.mesh -s example/fourlayer/FourLayer.source -p example/fourlayer/FourLayer.opt -m is  -t 32

time -p ./source/timos -f example/cube_5med/cube_5med.mesh -s example/cube_5med/cube_5med.source -p example/cube_5med/cube_5med.opt -m is  -t 32

time -p ./source/timos -f example/mouse/mouse.mesh -s example/mouse/mouse.source -p example/mouse/mouse.opt -m is  -t 32

time -p ./source/timos -f example/half_sphere/spherelens.mesh -s example/half_sphere/sphere.source -p example/half_sphere/freespace.opt -m is  -T 0.0001 200 -t 32
time -p ./source/timos -f example/half_sphere/spherelens.mesh -s example/half_sphere/sphere.source -p example/half_sphere/tissue.opt -m is  -T 0.0001 200 -t 32
