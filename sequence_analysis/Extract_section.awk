BEGIN {
	#print "Tile\tx\ty\tRead_alignment";
	split(pos_string, positions, " ")
}
{
	if (NR==1) {seq="Selection_sequence"; qual="Selection_quality"}
	else {
		seq = ""
		qual = ""
		for (i in positions) {
			p = positions[i]
			seq = seq substr($4, p, 1);
			qual = qual substr($5, p, 1);
		}
	}
	print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" seq "\t" qual	
}
