import sys
def sub(path, old, new):
    s = open(path, newline='').read()
    nl = '\r\n' if '\r\n' in s else '\n'
    old, new = old.replace('\n', nl), new.replace('\n', nl)
    assert s.count(old) >= 1, (path, old[:60])
    s = s.replace(old, new, 1)
    open(path, 'w', newline='').write(s)
d = sys.argv[1]; k = d + '/kepu.f'
sub(k, "              MFRAC(I) = MUF(I)/SIGMA\n           END DO\n        END IF\nc...  Define vectors & amplitudes (formulas differ depending on whether",
       "              MFRAC(I) = MUF(I)/SIGMA\n           END DO\n        ELSE\n           GM = MUF\n        END IF\nc...  Define vectors & amplitudes (formulas differ depending on whether")
sub(k, "           AMX = -SIGMA*SOM\n           AMY = SQRT(QQ*GM*(1.d0+EXC))*COM",
       "           AMX = -GM*SI*SOM\n           AMY = SQRT(QQ*GM*(1.d0+EXC))*SI*COM")
if len(sys.argv) > 2:
    sub(k, "          RAPQ = RAPQ*(Q/QT)*(UPQT/UPQ)*\n", "          RAPQ = RAPQ*\n")
print('ok', d)
