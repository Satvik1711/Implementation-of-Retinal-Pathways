[Transient Analysis]
{
   Npanes: 4
   Active Pane: 1
   {
      traces: 1 {524291,0,"V(vdiff2)"}
      X: ('m',0,0,0.003,0.03)
      Y[0]: ('m',0,-0.35,0.07,0.68)
      Y[1]: ('p',0,1e+308,9e-012,-1e+308)
      Volts: ('m',0,0,0,-0.35,0.07,0.68)
      Log: 0 0 0
      GridStyle: 1
      Text: "V" 1 (0.0958262831359278,1.00524585492228) ;time
      Text: "V" 1 (0.0285284280936455,-0.417823834196891) ;time
      Text: "V" 1 (0.0147658862876254,0.531502590673575) ;ON spikes
   },
   {
      traces: 1 {524293,0,"V(vlog)"}
      X: ('m',0,0,0.003,0.03)
      Y[0]: ('m',0,0.25,0.05,0.8)
      Y[1]: ('n',1,1e+308,1e-010,-1e+308)
      Volts: ('m',0,0,0,0.25,0.05,0.8)
      Log: 0 0 0
      GridStyle: 1
      Arrow: "V" 1 0 (0.0656514382402707,-0.642105263157895) (0.0606880992667795,-0.642105263157895)
      Circle: "V" 1 0 (0.000652173913043478,1.25157894736842) (0.0104849498327759,-0.421578947368421)
      Circle: "V" 1 0 (0.000652173913043478,1.25157894736842) (0.0104849498327759,-0.421578947368421)
      Text: "" 1 (0.0721376198533559,-0.636631578947369) ;Vph proportinal to ln(Iph)
      Text: "V" 1 (0.00501672240802676,0.58) ;ON spikes for positive slope of Vlog
      Text: "V" 1 (0.0219732441471572,0.727631578947369) ;OFF spikes for negative slope of Vlog
      Text: "V" 1 (0.0219732441471572,0.727631578947369) ;OFF spikes for negative slope of Vlog
      Text: "V" 1 (0.0219732441471572,0.727631578947369) ;OFF spikes for negative slope of Vlog
   },
   {
      traces: 1 {524298,0,"V(vdiff)"}
      X: ('m',0,0,0.003,0.03)
      Y[0]: ('m',0,-0.42,0.06,0.24)
      Y[1]: ('_',0,1e+308,9.9999999999999e-006,-1e+308)
      Volts: ('m',0,0,0,-0.42,0.06,0.24)
      Log: 0 0 1
      GridStyle: 1
      Arrow: "V" 1 0 (0.0719684151156232,9.40903260618677e-011) (0.067005076142132,9.40903260618677e-011)
      Circle: "V" 1 0 (0.0162541806020067,0.662539682539683) (0.0260869565217391,-1.35587301587302)
      Circle: "V" 1 0 (0.0162541806020067,0.662539682539683) (0.0260869565217391,-1.35587301587302)
      Text: "A" 1 (0.0743936830231246,2.65021082331517e-010) ;ln(Iph)
      Text: "V" 1 (0.0281605351170569,-0.00793650793650786) ;OFF spikes
   },
   {
      traces: 2 {589835,0,"V(vph)"} {34603010,1,"I(I1)"}
      X: ('m',0,0,0.003,0.03)
      Y[0]: ('m',0,0,0.04,0.48)
      Y[1]: ('n',0,0,1e-009,1.2e-008)
      Volts: ('m',0,0,0,0,0.04,0.48)
      Amps: ('n',0,0,0,0,1e-009,1.2e-008)
      Log: 0 0 0
      GridStyle: 1
      Text: "V" 1 (0.0131772575250836,0.266666666666667) ;Photocurrent
   }
}
