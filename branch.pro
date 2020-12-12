FUNCTION branch,sion

;returns array of barnching ratia for weak s path
;array element is defined to be 1d0 unless there is a branching
;branching ratios calculated from: Takahashi, K., Yokoi, K., At. Data Nucl. Data Tables 36, 375 (1987)

;branchings (defined by b+)
b=dblarr(n_elements(sion))
b[*]=1d0

Cu63=where(sion eq 'Cu63')
b[Cu63[0]]=0.3311d0
b[Cu63[1]]=1d0-b[Cu63[0]]

Br79=where(sion eq 'Br79')
b[Br79[0]]=0.0325d0
b[Br79[1]]=1d0-b[Br79[0]]



RETURN,b

END