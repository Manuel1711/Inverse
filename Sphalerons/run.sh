#!/bin/bash

set -e

PATH=$PATH:/Users/manuel/Documents/GitHub/Conducibility/Programmi/

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@             SMEARING EXTRACTION            @@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#


# This is a bash script that allows to solve inverse problems changing
#different parameters, such as the basis, the method and the simulation
#parameters. In the following everything is commented so you can
#easily re-use the code. This code takes inputs by inputfile.txt.
#So please change it.



optimization="0 -g"

#Choice of the method you want to use.
#If 0 -> Regular
#If 1 -> Modified
#method
boolM=1


#Choice of the basis function
#If 0 exp
#If 1 cosh
#basis
boolB=1
#declared="-D ABAB"
 
#Define Method 
if [ $boolM -eq 0 ]
then
    Method="-D BG"
fi

if [ $boolM -eq 1 ]
then
    Method="-D HLN"
fi


#Define Basis
if [ $boolB -eq 0 ]
then
    Basis="-D EXP"
fi

if [ $boolB -eq 1 ]
then
    Basis="-D COS"
fi

if [ $boolB -eq 2 ]
then
    Basis="-D COS_SPHAL"
fi


#Reading of the input file
reading="0"
parameters=()
par_name=( "Ls" "#Punti" "Nt" "Lside" "Rside" "Trash" "Sigma" "Apar")
filename=inputfile.txt
while read line; do
    if [[ "$line" == *"HIDDEN"* ]]; then
        # If the variable is already "1", It puts it equal to "0"
        if [[ "$reading" == "1" ]]; then
            reading="0"
        # Differently it puts it equal to "1"
        else
            reading="1"
        fi
        #Jumps to the next row
        continue
    fi
    if [[ "$reading" == "0" ]]; then
    # Jumps the rows that start with #
	if grep -Eq '^#|@' <<< "$line"; then
            continue
	fi
	# Extract only the numbers and it prints them
	number=$(echo "$line" | grep -Eo '[0-9]+(\.[0-9]+)?')
	parameters+=("$number")
    fi
done < "$filename"
 
echo "File di input:"
for (( i=0; i<${#par_name[@]}; i++ )); do
    echo "${par_name[i]}=${parameters[i]}"
done


#Compile with all the libraries
g++ -O$optimization -std=c++14 -o Smearing_Func Smearing_Func.C -I/Users/manuel/Desktop/gmpfrxx  -L/Users/manuel/Desktop/gmpfrxx -I/usr/local/include -L/usr/local/include  -lgmpfrxx -lmpfr -lgmpxx -lgmp -lm -lgsl  -lgslcblas -I/Users/manuel/Documents/GitHub/Inverse -L/Users/manuel/Documents/GitHub/Inverse $Method $Basis $N





:<<EOF 
if [ -e Output/Alpha/Ns60_Nt20/N_Cool_18/Output_ens.txt ]; then
	rm Output/Alpha/Ns60_Nt20/N_Cool_18/Output_ens.txt
fi
EOF

:<<EOF
valoriOmega=()
for i in {0..50..1}
do
    valore=$(bc <<< "scale=1; $i * 0.1")
    valoriOmega+=($valore)
done


for i in {1..50}
do



if [ -e Output/Alpha/Ns${parameters[0]}_Nt${parameters[2]}/Output_ens.txt ]; then
    rm Output/Alpha/Ns${parameters[0]}_Nt${parameters[2]}/Output_ens.txt
fi
EOF

echo Output/Alpha/Ns${parameters[0]}_Nt${parameters[2]}/Output_ens.txt

if [[ "${parameters[2]}" == "12" ]]; then
    sigma=(0.125 0.1333333333 0.1416666666666667 0.145833333333333 0.15 0.1583333333333333 0.1666666666666667 0.175 0.2083333333333333 0.25 0.2916666666666667 0.3333333333333333 0.375) 
elif [[ "${parameters[2]}" == "10" ]]; then
    sigma=(0.15 0.16 0.17 0.175 0.18 0.19 0.2 0.21 0.25 0.3 0.35 0.4 0.45) 
elif [[ "${parameters[2]}" == "14" ]]; then
    sigma=(0.1071428 0.1142857 0.12142857 0.125 0.12857 0.1357 0.142857 0.15 0.178571428 0.2142857 0.25 0.2857142857 0.32142857)
elif [[ "${parameters[2]}" == "16" ]]; then
    sigma=(0.09375 0.1 0.10625 0.109375 0.1125 0.11875 0.125 0.13125 0.15625 0.1875 0.21875 0.25 0.28125)  
elif [[ "${parameters[2]}" == "8" ]]; then
    sigma=(0.1875 0.2 0.2125 0.21875 0.225 0.2375 0.25 0.2625 0.3125 0.375 0.4375 0.5 0.625)
elif [[ "${parameters[2]}" == "1" ]]; then
    sigma=(1.5 1.6 1.7 1.75 1.8 1.9 2 2.1 2.5 3 3.5 4 4.5) 
fi

##Lunch the program making a cycle on the parameter of cooling
for i in {7..7}
#for ((i=6; i<=10; i+=2))
do
for j in {0..12}
do
./Smearing_Func S${parameters[0]} N${parameters[1]} B${parameters[2]} D${parameters[3]} U${parameters[4]} T${parameters[5]} E${sigma[j]} A${parameters[7]} L$i O0
done
done
 
:<<EOF
anche NBoot
./Smearing_Func L18 S60 N10 B20 O0 D35000.89 U1000.89 T1
./Smearing_Func L$i S48 N8 B16 O0 D30000.89 U50.89
./Smearing_Func L1 S36 N6 B12 O0 D10.92 U1.89
done
EOF
#A
