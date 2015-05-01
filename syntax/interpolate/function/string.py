def partition(expr, sep=(' ','\n', '\t', '\r')):
    b = 0
    result = []
    l = [len(s) for s in sep]
    n = len(sep)
    i = 0
    l1 = len(expr)
    while i<l1:
        for j in xrange(n):
	    if expr[i:i+l[j]] == sep[j]:
	        if b==i:
	            result.append(sep[j])
		else:
	            result.extend([expr[b:i], sep[j]])
		b = i+l[j]
		i = b-1
		break
	i += 1
    if l1-b>0:		
        result.append(expr[b:])		
    return result

def parseStr(x):
    if x.isalpha(): return x
    elif x.isdigit(): return int(x)
    elif x[-1] <> 'j':
        try:
            return float(x)
        except ValueError:
            return x
    else:
        try:
            return complex(x)
        except ValueError:
	    return x
    return x

