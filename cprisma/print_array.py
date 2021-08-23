def print_array(title,dic_type,ali_log):

  lent=len(title)
  list_f=[]
  for pri_v in dic_type.values():
      list_f.append(len(str(pri_v)))
  max_f=max(list_f)
  ax=0
  for pri_k, pri_v in dic_type.items():
      z=max_f-list_f[ax]
      if ax ==0:
          print(f"{title}| {pri_k} {pri_v}{' '*z} |")
          ali_log.write(f"{title}| {pri_k} {pri_v}{' '*z} |\n")
      else:
          print(f"{' '*lent}| {pri_k} {pri_v}{' '*z} |")
          ali_log.write(f"{' '*lent}| {pri_k} {pri_v}{' '*z} |\n")
      ax+=1
  ali_log.write(f"\n")
