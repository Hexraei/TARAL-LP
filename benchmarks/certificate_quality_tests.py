#!/usr/bin/env python3
"""Reporting-only quality and raw magnitude metric regression, optional GREENBEA."""
import argparse,json,subprocess,tempfile
from pathlib import Path
from simplex_certificate_check import verify

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--engine',required=True);ap.add_argument('--greenbea');a=ap.parse_args()
    checked=0
    with tempfile.TemporaryDirectory() as td:
        js=Path(td)/'out.json'
        def run(path,method):
            subprocess.run([a.engine,str(path),'--method',method,'--json',str(js)],check=True,capture_output=True,timeout=75)
            return json.loads(js.read_text())
        p=Path(td)/'tiny.mps';p.write_text('NAME TINY\nROWS\n N obj\n E zero\nCOLUMNS\n x obj 0 zero 1\nRHS\nBOUNDS\n FX b x 5e-10\nENDATA\n')
        for method in ['simplex','dual']:
            r=run(p,method);assert r['status']=='optimal' and r['certificate_quality']=='pass'
            assert r['row_violation_abs']==[5e-10] and r['row_term_magnitude']==[5e-10]
            assert r['row_violation_magnitude_scaled']==[1.] and verify(p,r)['pass_certificate'];checked+=1
            # A fake strict-quality label must not be trusted by the checker.
            r['certificate_quality']='unknown';assert not verify(p,r)['pass_certificate'];checked+=1
        p.write_text('NAME LIMIT\nROWS\n N obj\nCOLUMNS\n x obj -1\nENDATA\n')
        rc=subprocess.run([a.engine,str(p),'--time-limit','-1','--json',str(js)],capture_output=True).returncode;assert rc in (0,4)  # 4 = time_limit exit code on main
        assert json.loads(js.read_text())['certificate_quality']=='unknown';checked+=1
        if a.greenbea:
            for method in ['simplex','dual']:
                r=run(a.greenbea,method);v=verify(a.greenbea,r)
                assert r['status']=='optimal' and r['certificate_quality']=='fail'
                assert v['certificate_quality']=='fail' and not v['pass_certificate'] and not v['mismatched_fields']
                assert 1e-8 < r['primal_res'] < 1e-6
                for raw,den,scaled in zip(r['row_violation_abs'],r['row_term_magnitude'],r['row_violation_magnitude_scaled']):
                    assert abs(scaled-(raw/den if den else raw))<1e-15
                checked+=1
    print(f'PASS {checked} quality/reporting checks')
if __name__=='__main__':main()
