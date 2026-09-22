exec(open(__file__.replace('decode_rows.py','decode_probe.py'),encoding='utf8').read().split('for maxw')[0])
starts=np.frombuffer(b, '<u4',775,6580+n)
print('remaining:',b[6580+n+3100:].hex(' '));print('shorts',np.frombuffer(b[6580+n+3100:6580+n+3128],'<i2'))
for y in range(775):
 s=b[6580+starts[y]:6580+(starts[y+1] if y<774 else n)];vals=[]; i=0;sing=[]
 while i<len(s):
  c=s[i]; i+=1;lo=c&15;hi=c>>4
  if 1<=lo<=4 and 1<=hi<=4:
   for w in [lo,hi]:
    v=int.from_bytes(s[i:i+w],'little');i+=w
    vals.extend([(v>>(j*w)&((1<<w)-1))-((1<<(w-1))-1) for j in range(8)])
  else:
   vals.append(c-127)
   if c!=127:sing.append((i-1,c,len(vals)))
 if len(vals)!=800 or abs(sum(vals))>2:print('row',y,'n',len(vals),'sum',sum(vals),'singletons',sing,'hex',s.hex(' ') if len(vals)!=800 else '')
