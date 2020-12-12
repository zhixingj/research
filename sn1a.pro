function sn1a, one_a, z, z_max, b, num_iso, a

iso1A = DBLARR(num_iso,z_max)

f=tanh(b)
g=tanh(a-b)

for i = 0l,z_max-1 do begin
   for j = 0,num_iso-1 do begin
      iso1a(j,i)=(one_a(j))*z(i)*((tanh(a*z(i)-b))+f)/(g+f)
   endfor
endfor

RETURN, iso1a

END
