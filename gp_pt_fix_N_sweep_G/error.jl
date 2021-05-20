using DelimitedFiles

e1 = readdlm("error", comments=true,comment_char='#')

scale = 1000
relerr(A,B) = (A .- B)./B
mean(A) = sum(A)/length(A)

print("rspt ", scale*(relerr(e1[:,4].+e1[:,1], e1[:,2].+e1[:,1])), "\n");
print("enpt ", scale*(relerr(e1[:,5].+e1[:,1], e1[:,3].+e1[:,1])), "\n");
