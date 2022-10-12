[[ $* < 4 ]] && cat << AAA
## HELP: script mat_file component(e.g._Z for alpha, ZZ for beta) FIO1 FIO2 quality(optional, default 0)
## The script waits 3 seconds after printing the occupation numbers to give the possibility of cancelling the execution to change the parameters.
######
# Quality codes (from cubegen help in gaussian.com)
#  >> Number of points per side in the cube. A value of 0 selects the default value of 80**3 points distributed evenly over a rectangular grid generated automatically by the program (not necessarily a cube).
#  >> Positive values of npts similarly specify the number of points per side; e.g., 100 specified a grid of 1,000,000 (100**3) points.
#  >> The values -2, -3 and -4 correspond to the keywords Coarse, Medium and Fine and to values of 3 points/Bohr, 6 points/Bohr and 12 points/Bohr (respectively).
#  >> Negative values of npts <=-5 specify spacing of npts*10^-3 Angstroms between points in the grid.
######
AAA

quality_default=0
remove_ext() {                   ### Remove extension from a file
	echo "$1" | awk -F . 'BEGIN{sum=0};{for(i=1; i<NF;i++){sum+=length($i)}};END{ print substr($0,1,sum+NF-2)}'
}


file=${1:?mat_file component}
[ ! -f "${file}" ] && echo File not found && exit 1
comp=${2:?component}
comp="$(echo "${comp}" | tr [:lower:] [:upper:])"
[ "${#comp}" -lt 2 ] && comp="_${comp}"
[ "${#comp}" -gt 2 ] && echo -e "\n\n\t\tERROR: WRONG FIOs COMPONENT (${comp})\n" && exit 1
test="$(echo "${comp}" | tr -d '[X-Z]' | tr -d '[_]')"
[ "${#test}" -gt 0 ] && echo -e "\n\n\t\tERROR: WRONG FIOs COMPONENT (${comp})\n" && exit 1
alpha="-1"
[[ ${comp} =~ _ ]] && alpha="1"
fio1=${3:?FIO1}
fio2=${4:?FIO2}
case $fio1 in
    ''|*[!0-9]*) echo -e "\n\n\t\tERROR: WRONG FIO NUMBER ($fio1)\n" && exit 1;;
    *)  ;;
esac
case $fio2 in
    ''|*[!0-9]*) echo -e "\n\n\t\tERROR: WRONG FIO NUMBER ($fio2)\n" && exit 1;;
    *)  ;;
esac
compresor="$(which pbzip2)"
quality=${5:-${quality_default}}
echo -e "\n Mat file: ${file}, component: ${comp}, FIO1: ${fio1}, FIO2: ${fio2}, quality code:${quality}"
eval $(awk 	-v fio1="${fio1}" -v fio2="${fio2}" -v al="${alpha}"	'
	BEGIN									{
											a=0
											occ1=0
											occ2=0
										}
	/MO.* =.*   Vir/							{	nom=$2	}
	/ANALYSING '"${comp}"' FIELD COMPONENT OF HYPER\(POLARIZABILITY\)/ || /ANALYSING '"${comp}"' FIELD COMPONENT OF \(HYPER\)POLARIZABILITY/{a=1}
	a && /MO.*'"${fio1}"' occ =.*poptrans/					{	occ1=$5	}
	a && /MO.*'"${fio2}"' occ =.*poptrans/					{
											occ2=$5
											a=0
										}
	END									{
											if(nom<fio1)	{	print "OCC1=0"		}
											else		{	print "OCC1=" al*occ1	}
											if(nom<fio2)	{	print "OCC2=0"		}
											else		{	print "OCC2=" al*occ2	}
										}
								' "${file}")
[ "${OCC1}" = "0" -o "${OCC2}" = "0" ] && echo -e "\n\n\t\tERROR: the chosen FIOs do not exist (${fio1}, ${fio2})\n" && exit 1
echo -e " >> Substracting occupation numbers from ${file}: ${OCC1} and ${OCC2}"
[ "${alpha}" = "-1" ] && echo -e "    NOTE: for BETA the signs are changed."
fchkwo="$(remove_ext "${file}")_def${comp}"
fchk="${fchkwo}.fchk"
[ ! -f "${fchk}" ] && echo -e "\n\n\t\tERROR: fchk file (${fchk}) not found\n" && exit 1
cube1wo="$(remove_ext "${fchk}")-FIO${fio1}"
cube1="${cube1wo}.cube"
cube2wo="$(remove_ext "${fchk}")-FIO${fio2}"
cube2="${cube2wo}.cube"
echo -e "    YOU CAN CANCEL NOW THE SCRIPT IF THE FIOs PAIR IS NOT CORRECT"
sleep 3
echo -e " >> Writing cube 1 with quality code \"${quality}\": ${cube1} from ${fchk}"
cubegen 0 mo="${fio1}" "${fchk}" "${cube1}" "${quality}" h &> /dev/random
echo -e " >> Writing cube 2 copying quality: ${cube2} from ${fchk}"
cubegen 0 mo="${fio2}" "${fchk}" "${cube2}" -1 h "${cube1}" &> /dev/random

cubman_sq() {
	local in="${1:?}"
	local ou="${2:?}"
	cubman << EOF >& /dev/random
sq
${in}
y
${ou}
y
EOF
}

cube1sq="${cube1wo}-sq.cube"
cube2sq="${cube2wo}-sq.cube"
echo -e " >> Writing square of cube 1: ${cube1sq} from ${cube1}"
cubman_sq "${cube1}" "${cube1sq}"
echo -e " >> Writing square of cube 2: ${cube2sq} from ${cube2}"
cubman_sq "${cube2}" "${cube2sq}"

cubman_scale() {
	local in="${1:?}"
	local ou="${2:?}"
	local va="${3:?}"
	cubman << FIN >& /dev/random
sc
${in}
y
${ou}
y
${va}
FIN
}

cube1sc="${cube1wo}-sqxon.cube"
cube2sc="${cube2wo}-sqxon.cube"
echo -e " >> Multiplying by occupation number of cube 1: ${cube1sc} from ${cube1sq} x ${OCC1}"
cubman_scale "${cube1sq}" "${cube1sc}" "${OCC1}"
echo -e " >> Multiplying by occupation number of cube 2: ${cube2sc} from ${cube2sq} x ${OCC2}"
cubman_scale "${cube2sq}" "${cube2sc}" "${OCC2}"

cubman_add() {
	local in1="${1:?}"
	local in2="${2:?}"
	local out="${3:?}"
	cubman << ACABA >& /dev/random
a
${in1}
y
${in2}
y
${out}
y
ACABA
}

cubeadd="${fchkwo}_def_dens-${comp}-${fio1}_${fio2}.cube"
echo -e " >> Adding cube1 and cube2: ${cube2add} from ${cube1sc} + ${cube2sc}"
cubman_add "${cube1sc}" "${cube2sc}" "${cubeadd}"
echo -e " >> Compressing temporary files: ${cube1sc}, ${cube2sc}, ${cube1sq}, ${cube2sq}, ${cube1} and ${cube2}"
[ -f "${fchkwo}_def_dens-${comp}-${fio1}_${fio2}-cubes.tar.bz2" ] && rm "${fchkwo}_def_dens-${comp}-${fio1}_${fio2}-cubes.tar.bz2"
tar -I "$compresor" -cf "${fchkwo}_def_dens-${comp}-${fio1}_${fio2}-cubes.tar.bz2" "${cube1sc}" "${cube2sc}" "${cube1sq}" "${cube2sq}" "${cube1}" "${cube2}"
rm "${cube1sc}" "${cube2sc}" "${cube1sq}" "${cube2sq}" "${cube1}" "${cube2}"
echo -e "\n\t\t\t ... DONE ...\n"
