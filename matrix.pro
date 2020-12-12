FUNCTION matrix,uniq,sion,dion,sig,b

;CONSTRUCTS COEFFICIENT MATRIX
;ASSUMES THE USE OF COLUMN VECTORS FOR dN/dt = A*N
;IF DESIRE ROW VECTORS (ARRAYS) TAKE TRANSPOSE

a=DBLARR(n_elements(uniq),n_elements(uniq))

for i=0,n_elements(sion)-1 do begin
  m=where(uniq eq sion[i])
  a[m,m]=-1d0*sig[i]*b[i] + (-1d0*sig[i]*(1d0-b[i]))
endfor

for i=0,n_elements(dion)-1 do begin
  c=where(uniq eq dion[i])
  r=where(uniq eq sion[i])
  if c ne r then a[r,c]=sig[i]*b[i]
endfor

RETURN,a

END