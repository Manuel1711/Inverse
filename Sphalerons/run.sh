#!/bin/bash

set -e

PATH=$PATH:/Users/manuel/Documents/GitHub/Conducibility/Programmi/

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@             SMEARING EXTRACTION            @@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#


# This is a bash script that allows to solve inverse problems changing
#different parameters, such as the basis, the method and the simulation
#parameters. In the following everything is commented so you can
#easily re-use the code.



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
g++ -O$optimization -std=c++14 -o Smearing_Func Smearing_Func.C -I/Users/manuel/Desktop/gmpfrxx  -L/Users/manuel/Desktop/gmpfrxx -I/usr/local/include -L/usr/local/include  -lgmpfrxx -lmpfr -lgmpxx -lgmp -lm -lgsl  -lgslcblas -I/Users/manuel/Documents/GitHub/Inverse -L/Users/manuel/Documents/GitHub/Inverse  $Method $Basis $N





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

#Lunch the program making a cycle on the parameter of cooling
for i in {4..8}
do
./Smearing_Func S${parameters[0]} N${parameters[1]} B${parameters[2]} D${parameters[3]} U${parameters[4]} T${parameters[5]} E${parameters[6]} A${parameters[7]} L$i O0
done

:<<EOF
anche NBoot
./Smearing_Func L18 S60 N10 B20 O0 D35000.89 U1000.89 T1
./Smearing_Func L$i S48 N8 B16 O0 D30000.89 U50.89
./Smearing_Func L1 S36 N6 B12 O0 D10.92 U1.89
done
EOF
#A
