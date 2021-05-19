using DelimitedFiles

e1 = readdlm("error_l-0.50", comments=true,comment_char='#')
e2 = readdlm("error_l0.50", comments=true,comment_char='#')
e3 = readdlm("error_l-1.00", comments=true,comment_char='#')
e4 = readdlm("error_l1.00", comments=true,comment_char='#')

scale = 1000
relerr(A,B) = (A .- B)./A
mean(A) = sum(A)/length(A)

print("-0.50 rspt", scale*(relerr(e1[:,4], e1[:,2])), "\n");
print("-0.50 enpt", scale*(relerr(e1[:,5], e1[:,3])), "\n");

print(" 0.50 rspt", scale*(relerr(e2[:,4], e2[:,2])), "\n");
print(" 0.50 enpt", scale*(relerr(e2[:,5], e2[:,3])), "\n");

print("-1.00 rspt", scale*(relerr(e3[:,4], e3[:,2])), "\n");
print("-1.00 enpt", scale*(relerr(e3[:,5], e3[:,3])), "\n");

print(" 1.00 rspt", scale*(relerr(e4[:,4], e4[:,2])), "\n");
print(" 1.00 enpt", scale*(relerr(e4[:,5], e4[:,3])), "\n");
