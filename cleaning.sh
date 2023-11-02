#!/bin/bash


sed 's/[][",'"'"']//g' berry_averaged_points.txt | awk 'NF>1{print}' > berry_averaged_points.cleaned.txt
sed 's/[][",'"'"']//g' berry_detailed_points.txt | awk 'NF>1{print}' > berry_detailed_points.cleaned.txt
sed 's/[][",'"'"']//g' magnonic_surfaces.txt | awk 'NF>1{print}' > magnonic_surfaces.cleaned.txt
sed 's/[][",'"'"']//g' k_points.txt | awk 'NF>1{print}' > k_points.cleaned.txt
sed 's/[][",'"'"']//g' k_points_dynamically_refined.txt |  awk 'NF>1{print}' > k_points_dynamically_refined.cleaned.txt
