#!/bin/bash
PROGNAME="$(basename $0)"
export LC_ALL=C

yesno=("yes" "no")
quality_default=0
scaled_default=0
convert_default=0

fio1_default=1
fio2_default=2
compressor_default="pbzip2"

usage() {
cat << AAA

   HELP: Long  arguments: ${PROGNAME} --file mat_file --component "text" -1 "number" -2 "number" --quality --root-name "text" --scaled
         Short arguments: ${PROGNAME} -f mat_file -c -1 FIO1 -2 FIO2 -q -s
         Options:
             * -h|--help       :  this help
             * -f|--file       :  ".mat" file with EDOs/FIOs results. (mandatory)
             * -c|--component  :  component to be calculated. "text" corresponds to "x", "y" or "z" for alpha, "xx", ... for beta or "edo" for EDOs. (mandatory)
	     * -C|--compressor :  compressor to be used. (optional, default "${compressor_default}")
             * -1              :  "number" of the first  EDO/FIO. (optional, default "${fio1_default}")
             * -2              :  "number" of the second EDO/FIO. (optional, default "${fio2_default}")
             * -o|--convert    :  convert last cube file from MO to density format. (optional, default "${yesno[${convert_default}]}")
             * -q|--quality    :  cube files quality. (optional, default "${quality_default}")
             * -r|--root-name  :  change the root name of the fchk, following the "text", necessary to create the cube files.
             * -R|--rename     :  rename the last file with this name.
             * -s|--scaled     :  values to be scaled by the occupation number of the EDO/FIO. (optional, default "${yesno[${scaled_default}]}")
 
   The script waits 3 seconds after printing the occupation numbers to give the possibility of cancelling the execution to change the parameters.

   Quality codes (from cubegen help in gaussian.com)
    >> Number of points per side in the cube. A value of 0 selects the default value of 80**3 points distributed evenly over a rectangular grid generated automatically by the program (not necessarily a cube).
    >> Positive values of npts similarly specify the number of points per side; e.g., 100 specified a grid of 1,000,000 (100**3) points.
    >> The values -2, -3 and -4 correspond to the keywords Coarse, Medium and Fine and to values of 3 points/Bohr, 6 points/Bohr and 12 points/Bohr (respectively).
    >> In order to specify a reference cube file, type the filename after the argument "-q/--quality".
    >> Negative values of npts <=-5 specify spacing of npts*10^-3 Angstroms between points in the grid.

AAA
}

remove_ext() {                   ### Remove extension from a file name
	echo "$1" | awk -F . 'BEGIN{sum=0};{for(i=1; i<NF;i++){sum+=length($i)}};END{ print substr($0,1,sum+NF-2)}'
}

extract_ext() {                   ### Extract extension from a file name
	echo ${1} | awk -F. ' { print $NF }'
}

extract_2ext() {                   ### Extract double extension from a file name
	echo "${1}" | awk -F. ' { print $(NF-1) "." $NF }'
}

def_var() {
	local name ifnotfound varname
	name="${1}"
	ifnotfound="${2}"  # =1 ... exit
	varname="${3}"
	[[ "${#varname}" -eq 0 ]] && varname="${name}"
	if which "${name}" > /dev/null 2>&1
		then
			eval ${varname}="$(which ${name})"
		else
			if [ "${ifnotfound}" != 0 ]
				then
					echo -e "\n\n\t\tERROR: ${name} WAS NOT FOUND\n\t\tType \"-h\" or \"--help\" to output a usage message and exit.\n" && exit "${ifnotfound}"
				else
					eval ${varname}=""
			fi
	fi
}

reverse() {
        local def=$1
        local var=$2
        local val
        eval val=\${${def}}
        [ "${val}" = "0" ] && eval ${var}="1" || eval ${var}="0"
}

scaled="${scaled_default}"
convert="${convert_default}"
rootname=""
rename=""
GETOPT_ARGS=$(getopt -o hf:1:2:oq:sc:C:r:R: -l "compressor:","component:","help","file:","quality:","scaled","convert","root-name:","rename:" -n "${PROGNAME}" -- "$@")
[[ $? -ne 0 ]] && echo -e "\n\n\t\tERROR: UNEXPECTED OR UNKNOWN ARGUMENT\n\t\tType \"-h\" or \"--help\" to output a usage message and exit.\n" && exit 1
eval set -- "$GETOPT_ARGS"

while :; do
	case "$1" in
		-c|--component)
			shift
			comp="$1"
			shift
			;;
		-C|--compressor)
			shift
			compressor="$1"
			shift
			;;
		-h|--help)
			usage
			exit 0
			;;
		-f|--file)
			shift
			file="$1"
			shift
			;;
		-1       )
			shift
			fio1="$1"
			shift
			;;
		-2       )
			shift
			fio2="$1"
			shift
			;;
		-o|--convert)
			reverse convert_default convert
			shift
			;;
		-q|--quality)
			shift
			quality="$1"
			shift
			;;
		-r|--root-name)
			shift
			rootname="$1"
			shift
			;;
		-R|--rename)
			shift
			rename="$1"
			shift
			;;
		-s|--scaled)
			reverse scaled_default scaled
			shift
			;;
		--)
			shift
			break
			;;
	esac
done

def_var cubegen 1
def_var cubman  1
compressor="${compressor:-${compressor_default}}"
def_var "${compressor}" 0 compressor
[[ ${#compressor} -gt 0 ]] && compressor="-I ${compressor}"
[[ "${#file}" -lt 1 ]] && echo -e "\n\n\t\tERROR: you must specify a \"file\"\n\t\tType \"-h\" or \"--help\" to output a usage message and exit.\n" && exit 1
[ ! -f "${file}" ] && echo -e "\n\n\t\tERROR: FILE NOT FOUND\n\t\tType \"-h\" or \"--help\" to output a usage message and exit.\n" && exit 1
[[ ${#comp} -lt 1 ]] && echo -e "\n\n\t\tERROR: COMPONENT NOT SPECIFIED\n\t\tType \"-h\" or \"--help\" to output a usage message and exit.\n" && exit 1
comp="$(echo "${comp}" | tr [:lower:] [:upper:])"
if [ "${comp}" != "EDO" ]
	then
		[ "${#comp}" -gt 2 ] && echo -e "\n\n\t\tERROR: WRONG FIOs COMPONENT (${comp})\n\t\tType \"-h\" or \"--help\" to output a usage message and exit.\n" && exit 1
		test="$(echo "${comp}" | tr -d '[X-Z]' | tr -d '[_]')"
		[ "${#test}" -gt 0 ] && echo -e "\n\n\t\tERROR: WRONG FIOs COMPONENT (${comp})\n\t\tType \"-h\" or \"--help\" to output a usage message and exit.\n" && exit 1
		text="FIO"
	else
		text="EDO"
fi
scaled="${scaled:-${scaled_default}}"
convert="${convert:-${convert_default}}"
quality="${quality:-${quality_default}}"
fio1="${fio1:-${fio1_default}}"
fio2="${fio2:-${fio2_default}}"

alpha="1"
[[ "${#comp}" == 2 ]] && alpha="-1"
case ${fio1} in
    ''|*[!0-9]*) echo -e "\n\n\t\tERROR: WRONG ${text} NUMBER ($fio1)\n" && exit 1;;
    *)  ;;
esac
case ${fio2} in
    ''|*[!0-9]*) echo -e "\n\n\t\tERROR: WRONG ${text} NUMBER ($fio2)\n" && exit 1;;
    *)  ;;
esac

unset ref_cubefile

is_integer() {
    return $(test "$@" -eq "$@" > /dev/null 2>&1);
}

if is_integer "${quality}"
	then
		[[ "${quality}" -eq "-1" ]] && echo -e "\n\n\t\tERROR: Wrong quality code. In order to specify a cube file name, write it directly\n" && exit 1
	else
		[ ! -f "${quality}" ] && echo -e "\n\n\t\tERROR: Reference cube file (\""${quality}"\") NOT FOUND\n" && exit
		ref_cubefile="${quality}"
		quality="-1"
		if [ "$(extract_ext $ref_cubefile)" != "cube" ]
			then
				ref_cube_compress="echo"
				ref_cube_uncompress="echo"
				if [ "$(extract_2ext $ref_cubefile)" = "cube.gz"  ] 
					then
						ref_cube_compress="gzip"
						ref_cube_uncompress="gunzip"
				elif [ "$(extract_2ext $ref_cubefile)" = "cube.bz2" ]
					then
						ref_cube_compress="bzip2"
						ref_cube_uncompress="bunzip2"
				else
					echo -e "\n\n\t\tERROR: Unknown compressed cube file extension\n" && exit 1
				fi
				trap "${ref_cube_compress} $(remove_ext ${ref_cubefile}) &> /dev/random" EXIT
				"${ref_cube_uncompress}" "${ref_cubefile}" &> /dev/random
				ref_cubefile="$(remove_ext "${ref_cubefile}")"
		fi
fi

if [ "${text}" = "EDO" ] 
	then
		eval $(awk '/ANALYSING .* FIELD COMPONENT OF ELECTRON DEFORMATION ORBITALS/{print "comp=" $3}' "${file}")
		comp="$(echo "${comp}" | tr [:upper:] [:lower:])"
		fchkwo="$(remove_ext "${file}")_${text}-${comp}"
	else
		fchkwo="$(remove_ext "${file}")_${text}${comp}"
fi
[[ "${#rootname}" -gt 0 ]] && fchkwo="${rootname}"
fchk="${fchkwo}.fchk"
[ ! -f "${fchk}" ] && echo -e "\n\n\t\tERROR: fchk file (${fchk}) not found\n" && exit 1
echo -e "\n Mat file: ${file}, component: ${comp}, ${text}1: ${fio1}, ${text}2: ${fio2}, quality code:${quality}"
eval $(awk 	-v fio1="${fio1}" -v fio2="${fio2}" -v al="${alpha}" -v sca="${scaled}"	'
	function abs(v)   {return v < 0 ? -v : v}
	BEGIN									{
											a=0
											occ1=0
											occ2=0
										}
	/MO *[0-9]*.* = .*/							{	nom=$2	}
	/ANALYSING '"${comp}"' FIELD COMPONENT OF HYPERPOLARIZABILITY/ || /ANALYSING '"${comp}"' FIELD COMPONENT OF POLARIZABILITY/ || /ANALYSING .* FIELD COMPONENT OF ELECTRON DEFORMATION ORBITALS/{a=1}
	a && /'"${text}"' *'"${fio1}"' occ =.*poptrans/					{	occ1=$5	}
	a && /'"${text}"' *'"${fio2}"' occ =.*poptrans/					{
											occ2=$5
											a=0
										}
	END									{
											if(nom<fio1)	{	print "OCC1=0"			}
											else if(sca!=0)	{	print "OCC1=" al*occ1/abs(occ1)	}
											else		{	print "OCC1=" al*occ1		}
											if(nom<fio2)	{	print "OCC2=0"			}
											else if(sca!=0)	{	print "OCC2=" al*occ2/abs(occ2)	}
											else		{	print "OCC2=" al*occ2		}
										}
								' "${file}")
[ "${OCC1}" = "0" -o "${OCC2}" = "0" ] && echo -e "\n\n\t\tERROR: the chosen EDOs/FIOs do not exist (${fio1}, ${fio2})\n" && exit 1
echo -e " >> Occupation numbers from ${file}: ${OCC1} and ${OCC2}"
[ "${alpha}" = "-1" ] && echo -e "    NOTE: for BETA the signs are changed."
cube1wo="$(remove_ext "${fchk}")-${fio1}"
cube1="${cube1wo}.cube"
cube2wo="$(remove_ext "${fchk}")-${fio2}"
cube2="${cube2wo}.cube"
echo -e "    YOU CAN CANCEL NOW THE SCRIPT IF THE ${text}s PAIR IS NOT CORRECT"
sleep 3
echo -en " >> Writing cube 1 with quality code \"${quality}\": ${cube1} from ${fchk}"
[ ! -z "$ref_cubefile" ] && echo ". Using grid from \"${ref_cubefile}\"" || echo -e "\n"
"${cubegen}" 0 mo="${fio1}" "${fchk}" "${cube1}" "${quality}" h "${ref_cubefile}" &> /dev/random
[ ! -z "$ref_cubefile" ] && "${ref_cube_compress}" "${ref_cubefile}"
echo -e " >> Writing cube 2 copying quality: ${cube2} from ${fchk}"
"${cubegen}" 0 mo="${fio2}" "${fchk}" "${cube2}" -1 h "${cube1}" &> /dev/random

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

convert_cube() {
	local incube="${1:?}"
	local tmpcube="tmp$(date +%g%m%d%H%M%S).cube"
	awk '
		BEGIN 	{ 	pos=0 }
		NR==1 	{ 	print }
		NR==2 	{ 	print "Electron density from Total SCF Density" }
		NR==3 	{
				nat=$1
				val1=$2
				val2=$3
				val3=$4
				if(nat<0){
						pos=0
						nat=-nat
					 }
				printf "%5s%12s%12s%12s\n",nat,val1,val2,val3
			}
		NR>3&&NR<=nat+6
		pos&&NR==nat+7
		NR>nat+7
	    ' "${incube}" > "${tmpcube}"
	mv "${tmpcube}" "${incube}"
}

cubeadd="${fchkwo}_def_dens-${comp}-${fio1}_${fio2}.cube"
echo -e " >> Adding cube1 and cube2: ${cubeadd} from ${cube1sc} + ${cube2sc}"
cubman_add "${cube1sc}" "${cube2sc}" "${cubeadd}"
if [[ "${convert}" -eq 0 ]] 
	then
		echo -e " >> Converting cube type from MO to density file"
		convert_cube "${cubeadd}"
fi
if [[ "${#rename}" -gt 0 ]] 
	then
		echo -e " >> Renaming file from ${cubeadd} to ${rename}.cube"
		mv "${cubeadd}" "${rename}.cube"
fi
echo -e " >> Compressing temporary files: ${cube1sc}, ${cube2sc}, ${cube1sq}, ${cube2sq}, ${cube1} and ${cube2}"
[ -f "${fchkwo}_def_dens-${comp}-${fio1}_${fio2}-cubes.tar.bz2" ] && rm "${fchkwo}_def_dens-${comp}-${fio1}_${fio2}-cubes.tar.bz2"
tar ${compressor} -cf "${fchkwo}_def_dens-${comp}-${fio1}_${fio2}-cubes.tar.bz2" "${cube1sc}" "${cube2sc}" "${cube1sq}" "${cube2sq}" "${cube1}" "${cube2}"
rm "${cube1sc}" "${cube2sc}" "${cube1sq}" "${cube2sq}" "${cube1}" "${cube2}"
echo -e "\n\t\t\t ... DONE ...\n"
