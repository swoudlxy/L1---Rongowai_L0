% This function computes the distance from Tx to Sx to Rx
% based on their coordinates given in ECEF

function dis_abc = pdis(a,b,c)

dis1 = norm(c-a);   %tx2sx
dis2 = norm(b-c);   %sx2rx
dis_abc = dis1+dis2;