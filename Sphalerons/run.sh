#!/bin/bash

set -e

PATH=$PATH:/Users/manuel/Documents/GitHub/Conducibility/Programmi/

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
#@@@@@@@@@@@@             SMEARING EXTRACTION            @@@@@@@@@@@@#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#


#@@@ COMPLILAZIONE Spectral_gsl.C @@@#

optimization="0 -g"


#method
boolM=1
#basis
boolB=1
#declared="-D ABAB"
 

if [ $boolM -eq 0 ]
then
    Method="-D BG"
fi

if [ $boolM -eq 1 ]
then
    Method="-D HLN"
fi


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

reading="0"
parameters=()
par_name=( "Ls" "#Punti" "Nt" "Lside" "Rside" "Trash" "Sigma" "Apar")
filename=inputfile.txt
while read line; do
    if [[ "$line" == *"HIDDEN"* ]]; then
        # Se la variabile è già impostata su "1", la reimposta su "0"
        if [[ "$reading" == "1" ]]; then
            reading="0"
        # Altrimenti, imposta la variabile su "1"
        else
            reading="1"
        fi
        # Salta alla prossima riga
        continue
    fi
    if [[ "$reading" == "0" ]]; then
    # Salta le righe che iniziano con #
	if grep -Eq '^#|@' <<< "$line"; then
            continue
	fi
	# Estrae solo i numeri da ogni riga e li stampa
	number=$(echo "$line" | grep -Eo '[0-9]+(\.[0-9]+)?')
	parameters+=("$number")
    fi
done < "$filename"
 
echo "File di input:"
for (( i=0; i<${#par_name[@]}; i++ )); do
    echo "${par_name[i]}=${parameters[i]}"
done

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
