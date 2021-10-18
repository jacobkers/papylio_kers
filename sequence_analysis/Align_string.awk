BEGIN {
	print "Tile\tx\ty\tRead_aligned\tQuality_aligned";
}
$3==refname {
	#print $6;
	seq = "";
	qual = "";
	start = $4;
	for (i = 1; i <= start-1; ++i) {seq=seq "-"; qual=qual "-"};
	split($6, len, "M|I|D|N|S|H|P|=|X");
	split($6, op, "[0-9]*");
	#print $10

	curpos = 1;
	#print curpos;
	#for (bi in b) {print bi};
	for (i = 1; i <= length(op)-1; ++i) {
		if (op[i+1]=="M" || op[i+1]=="=" || op[i+1]=="X") {
			seq = seq substr($10, curpos, len[i]);
			qual = qual substr($11, curpos, len[i]);
			curpos = curpos + len[i];
		}
		else if (op[i+1]=="I") {curpos = curpos + len[i];}
		else if (op[i+1]=="S") {
			curpos = curpos + len[i];
			if (i>1) {for (j = 1; j <= len[i]; ++j) {seq=seq "-"; qual=qual "-"};};
		}
		else if (op[i+1]=="D" || op[i+1]=="N") {
			# if (i==1 && op[i+1]=="S") {len[i]=len[i]-curpos+1};
			for (j = 1; j <= len[i]; ++j) {seq=seq "-"; qual=qual "-"};
			# curpos = curpos + len[i];
		};

		#print curpos;
	};
	#print refseq
	reflen = length(refseq);
	seq = substr(seq,1,reflen);
	qual = substr(qual,1,reflen);

	split($1,a,":");
	print a[5] "\t" a[6] "\t" a[7] "\t" seq "\t" qual;
	#print a[5] "\t" a[6] "\t" a[7] "\t" refseq;
}
