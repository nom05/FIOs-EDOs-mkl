file="${1?mat file}"
[ ! -f "${file}" ] && echo File not found && exit 1
awk 	'
		function max(v,w)	{return v < w ?  w : v}
		function abs(v)		{return v < 0 ? -v : v}
		BEGIN			{
						n=0
						m=0
						c=0
						d=0
					}
		/HyperPolar =/		{
						for(i=3;i<=NF-2;i++)	{
										n++
										a[n]=$i
									}
					}
		/RE-COMPUTED BETA TENS./{
						for(i=4;i<=NF;i++)	{
										m++
										b[m]=$i
									}
					}
		END			{
						printf "a=\t"
						for(i=1;i<=n;i++)	{
										printf "%.4f ",a[i]
									}
						printf "\nb=\t"
						for(i=1;i<=n;i++)	{
										printf "%.4f ",b[i]
										c=max(c,abs(a[i]-b[i]))
										d=max(d,abs(abs(a[i])-abs(b[i])))
									}
						print "\n"
						printf "%.4f\t%.4f\n",c,d
					}
	'	 "${file}"
