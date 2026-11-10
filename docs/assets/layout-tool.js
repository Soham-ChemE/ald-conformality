/* Layout-to-recipe tool for the live page (section #layout-tool, markup in docs/partials/layout-tool.html).
   Classic script, no dependencies. Exposes window.LayoutTool = { compute, parseCsv, toCsv, presets, init }.
   compute() and parseCsv() are pure (no DOM) so tests/test_layout_tool.py can drive them headlessly.

   Physics, all from this repository:
   - supply rule   t_pulse = 1.2 (1 + A_wall/A_wafer) N0 / (U c0), A_wall/A_wafer = holeDensity pi d AR d per zone,
                   c0 = x p / (kB T), U = Q_sccm 1.6667e-8 (101325/p) (T/273.15) / (pi 0.15^2)   (same arithmetic as #try)
   - non-uniform layouts (paper/assumptions.md, Night 3 part 1): the wafer-average wall ratio under-predicts; the local ratio of
     the most demanding zone governs; the two chamber runs saturated at 0.82x (dense centre, 11.85/14.39) and 0.87x (AR ramp,
     12.50/14.36) of the local rule, 0.98x on the uniform wafer (14.05/14.39). That band is a calibration on two runs, not a model.
   - open field next to a dense array saturates in 0.40-0.50 s against 0.35 s on a blanket wafer (15-30 % slower). */
(function(){
  'use strict';
  var KB=1.380649e-23, WAFER_R_MM=150;
  var MARGIN={nonuniform:[0.82,0.87], uniform:[0.98,0.98]};
  var VALID={arMin:20, arMax:100, tmaMin:0.5, tmaMax:10};
  var FLAT_NEXT_TO_ARRAY=[0.40,0.50], FLAT_BLANKET=0.35;

  // ---------- presets ----------
  var presets={
    dense:{name:'Dense centre', zones:[
      {rIn:0,  rOut:75,  d:100, AR:100, density:1.78e13},
      {rIn:75, rOut:100, d:100, AR:100, density:8.48e12}, // linear ramp 1.78e13 -> 0 over 75-100 mm, replaced by its area-weighted mean (0.476x)
      {rIn:100,rOut:150, d:100, AR:0,   density:0}]},
    uniform:{name:'Uniform AR 100', zones:[{rIn:0, rOut:150, d:100, AR:100, density:1.78e13}]},
    ramp:{name:'AR 100 to 60 ramp', zones:[100,90,80,70,60].map(function(ar,i){return {rIn:30*i, rOut:30*(i+1), d:100, AR:ar, density:1.78e13};})}
  };
  var chamberDefault={x:1, Q:300, p:300, T:300, N0:4e18};

  // ---------- pure arithmetic ----------
  function supply(ch){
    var x=ch.x/100, p=ch.p, T=ch.T+273.15;
    var c0=x*p/(KB*T), U=ch.Q*1.6667e-8*(101325/p)*(T/273.15)/(Math.PI*.15*.15);
    return {c0:c0, U:U, ts0:ch.N0/(U*c0)};
  }
  function wallRatio(z){ var d=z.d*1e-9; return (z.density||0)*Math.PI*d*z.AR*d; }

  function compute(zones, chamber){
    var ch={}; for(var k in chamberDefault) ch[k]=(chamber&&isFinite(chamber[k]))?+chamber[k]:chamberDefault[k];
    var waferR=(chamber&&chamber.waferR)||WAFER_R_MM, waferArea=Math.PI*waferR*waferR;
    var S=supply(ch), warnings=[], rows=[], covered=0, wallArea=0;
    var zs=(zones||[]).map(function(z){return {rIn:+z.rIn, rOut:+z.rOut, d:+z.d, AR:+z.AR, density:+z.density||0};})
      .filter(function(z){return isFinite(z.rIn)&&isFinite(z.rOut)&&z.rOut>z.rIn&&z.rIn>=0;})
      .sort(function(a,b){return a.rIn-b.rIn;});
    for(var i=0;i<zs.length;i++){
      var z=zs[i], bare=!(z.density>0&&z.AR>0&&z.d>0), wall=bare?0:wallRatio(z);
      var ro=Math.min(z.rOut,waferR), area=Math.PI*(ro*ro-Math.min(z.rIn,waferR)*Math.min(z.rIn,waferR));
      if(i>0&&z.rIn<zs[i-1].rOut-1e-9) warnings.push('zones '+i+' and '+(i+1)+' overlap; the average counts the overlap twice');
      if(z.rOut>waferR+1e-9) warnings.push('zone '+(i+1)+' extends past the '+waferR+' mm wafer edge; clipped in the average');
      if(!bare&&(z.AR<VALID.arMin||z.AR>VALID.arMax)) warnings.push('zone '+(i+1)+': AR '+z.AR+' is outside the validated AR 20 to 100');
      covered+=area; wallArea+=wall*area;
      rows.push({rIn:z.rIn, rOut:z.rOut, d:z.d, AR:z.AR, density:z.density, bare:bare, wall:wall, area:area*1e-6, areaFrac:area/waferArea,
        local:1.2*(1+wall)*S.ts0});
    }
    if(covered<waferArea*(1-1e-6)) warnings.push(((1-covered/waferArea)*100).toFixed(0)+' % of the wafer is not in any zone; it is counted as bare in the average');
    if(ch.x<VALID.tmaMin||ch.x>VALID.tmaMax) warnings.push('TMA '+ch.x+' % is outside the validated 0.5 to 10 %');
    var holeRows=rows.filter(function(r){return !r.bare;});
    var gov=-1, govLocal=-Infinity;
    rows.forEach(function(r,i){ if(!r.bare&&r.local>govLocal){govLocal=r.local; gov=i;} });
    var avgWall=wallArea/waferArea, avgRule=1.2*(1+avgWall)*S.ts0;
    var walls=holeRows.map(function(r){return r.wall;});
    var uniform=holeRows.length>0&&covered>=waferArea*(1-1e-3)&&holeRows.length===rows.length&&
      (Math.max.apply(null,walls)-Math.min.apply(null,walls))<=1e-2*Math.max.apply(null,walls);
    var m=uniform?MARGIN.uniform:MARGIN.nonuniform;
    rows.forEach(function(r,i){
      r.governing=(i===gov);
      if(r.bare){ r.expectedLo=FLAT_NEXT_TO_ARRAY[0]; r.expectedHi=FLAT_NEXT_TO_ARRAY[1]; r.expectedNote='flat surface next to an array, from the dense-centre run; a blanket wafer takes '+FLAT_BLANKET.toFixed(2)+' s'; }
      else if(r.governing){ r.expectedLo=r.local*m[0]; r.expectedHi=r.local*m[1]; r.expectedNote=uniform?'uniform wafer: 0.98x the rule in the reference run':'0.82 to 0.87x the local rule in the two non-uniform runs'; }
      else { r.expectedLo=null; r.expectedHi=null; r.expectedNote='saturates before the recipe pulse ends; the local rule is not calibrated next to a denser neighbour'; }
    });
    var ars=holeRows.map(function(r){return r.AR;}), arMin=ars.length?Math.min.apply(null,ars):null, arMax=ars.length?Math.max.apply(null,ars):null;
    var recipe=gov>=0?{local:govLocal, lo:govLocal*m[0], hi:govLocal*m[1], zone:gov}:{local:null, lo:null, hi:null, zone:-1};
    return {chamber:ch, supply:S, zones:rows, uniform:uniform, margin:m, governing:gov, recipe:recipe,
      average:{wall:avgWall, rule:avgRule}, underDose:gov>=0?govLocal/avgRule:null,
      arSpread:ars.length?{min:arMin, max:arMax, ratio:arMax/arMin, timeSpreadPct:gov>=0?100*(1-(1+wallRatio({d:rows[gov].d,AR:arMin,density:rows[gov].density}))/(1+rows[gov].wall)):0}:null,
      warnings:warnings, valid:VALID};
  }

  // ---------- CSV ----------
  // lines "r_in_mm, r_out_mm, d_nm, AR, density_per_m2"; blank lines, '#' comments and a header line are skipped
  function parseCsv(text){
    var zones=[], errors=[], first=true;
    String(text||'').split(/\r?\n/).forEach(function(line,n){
      var s=line.trim(); if(!s||s[0]==='#') return;
      var f=s.split(/[,;\t]+/).map(function(v){return v.trim();}), header=first&&/^[a-z_]/i.test(f[0]); first=false;
      if(header) return;
      if(f.length<5){ errors.push('line '+(n+1)+': expected 5 fields'); return; }
      var v=f.slice(0,5).map(Number);
      if(v.some(function(q){return !isFinite(q);})){ errors.push('line '+(n+1)+': not a number'); return; }
      if(v[1]<=v[0]){ errors.push('line '+(n+1)+': r_out must exceed r_in'); return; }
      zones.push({rIn:v[0], rOut:v[1], d:v[2], AR:v[3], density:v[4]});
    });
    return {zones:zones, errors:errors};
  }
  function toCsv(zones){
    return ['r_in_mm, r_out_mm, d_nm, AR, density_per_m2'].concat((zones||[]).map(function(z){
      return [z.rIn, z.rOut, z.d, z.AR, (+z.density>0?(+z.density).toExponential(2).replace('e+','e'):'0')].join(', ');})).join('\n');
  }

  // ---------- rendering ----------
  var VIRIDIS=['#440154','#472d7b','#3b528b','#2c728e','#21918c','#27ad81','#5ec962','#aadc32','#fde725'];
  function hex2(h){return [parseInt(h.slice(1,3),16),parseInt(h.slice(3,5),16),parseInt(h.slice(5,7),16)];}
  function viridis(t){ t=Math.max(0,Math.min(1,isFinite(t)?t:0)); var k=t*(VIRIDIS.length-1), i=Math.floor(k), f=k-i; if(i>=VIRIDIS.length-1) return VIRIDIS[VIRIDIS.length-1];
    var a=hex2(VIRIDIS[i]), b=hex2(VIRIDIS[i+1]); return 'rgb('+a.map(function(c,j){return Math.round(c+(b[j]-c)*f);}).join(',')+')'; }
  function esc(s){return String(s).replace(/[&<>"]/g,function(c){return {'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c];});}
  function f2(v){return v==null||!isFinite(v)?'—':(+v).toFixed(2);}
  function f1(v){return v==null||!isFinite(v)?'—':(+v).toFixed(1);}
  function fden(v){ if(!(v>0)) return '0'; var e=Math.floor(Math.log10(v)), m=v/Math.pow(10,e); return m.toFixed(2)+' × 10<sup>'+e+'</sup>'; }
  function rangeTxt(lo,hi){ return lo==null?'—':(Math.abs(hi-lo)<5e-3?f2(lo):f1(lo)+'–'+f1(hi)); }

  function svgMap(res){
    var W=480, H=360, C=172, R=150, sc=R/WAFER_R_MM, maxT=res.recipe.local||1, parts=[], LX=352;
    parts.push('<circle cx="'+C+'" cy="'+C+'" r="'+R+'" fill="var(--line)" stroke="var(--muted)" stroke-width="1.5"/>');
    var cols=res.zones.map(function(z){return z.bare?'var(--bg2)':viridis(z.local/maxT);});
    res.zones.slice().reverse().forEach(function(z,ri){
      var i=res.zones.length-1-ri, ri0=Math.min(z.rIn,WAFER_R_MM)*sc, ro=Math.min(z.rOut,WAFER_R_MM)*sc, mid=(ri0+ro)/2, w=ro-ri0; if(w<=0) return;
      var tip='zone '+(i+1)+': '+z.rIn+'\u2013'+z.rOut+' mm'+(z.bare?', bare':', d '+z.d+' nm, AR '+z.AR+', '+(+z.density).toExponential(2)+' /m\u00b2, wall ratio '+f1(z.wall)+', local rule '+f2(z.local)+' s');
      parts.push('<circle cx="'+C+'" cy="'+C+'" r="'+mid+'" fill="none" stroke="'+cols[i]+'" stroke-width="'+w+'"'+(z.governing?' data-gov="1"':'')+'><title>'+esc(tip)+'</title></circle>');
    });
    parts.push('<circle cx="'+C+'" cy="'+C+'" r="'+R+'" fill="none" stroke="var(--muted)" stroke-width="1.5" pointer-events="none"/>');
    parts.push('<rect x="'+(C-10)+'" y="'+(C+R-3)+'" width="20" height="3" fill="var(--bg)" pointer-events="none"/>');
    // ring boundaries and zone numbers along the horizontal radius (alternating above and below the axis so narrow rings do not collide)
    res.zones.forEach(function(z,i){
      var ri0=Math.min(z.rIn,WAFER_R_MM)*sc, ro=Math.min(z.rOut,WAFER_R_MM)*sc, mid=C+(ri0+ro)/2; if(ro<=ri0) return;
      if(z.rIn>0) parts.push('<circle cx="'+C+'" cy="'+C+'" r="'+ri0+'" fill="none" stroke="var(--bg)" stroke-width="1" pointer-events="none"/>');
      parts.push('<text x="'+mid+'" y="'+(C+4+(i%2?14:-14))+'" text-anchor="middle" font-size="11" font-weight="700" fill="var(--ink)" stroke="var(--card)" stroke-width="3" paint-order="stroke" pointer-events="none">'+(i+1)+'</text>');
    });
    [50,100,150].forEach(function(r){ parts.push('<text x="'+(C+3)+'" y="'+(C-r*sc+11)+'" font-size="9" fill="var(--muted)" stroke="var(--card)" stroke-width="3" paint-order="stroke" pointer-events="none">'+r+' mm</text>'); });
    // zone list beside the wafer: swatch, radii, local rule
    parts.push('<text x="'+LX+'" y="22" font-size="10.5" font-weight="700" fill="var(--muted)">zone · r, mm · local rule</text>');
    res.zones.forEach(function(z,i){ var y=42+i*22; if(y>H-60) return;
      parts.push('<rect x="'+LX+'" y="'+(y-10)+'" width="12" height="12" rx="3" fill="'+cols[i]+'" stroke="var(--line)"/>');
      parts.push('<text x="'+(LX+18)+'" y="'+y+'" font-size="11" fill="'+(z.governing?'var(--acc2)':'var(--ink)')+'"'+(z.governing?' font-weight="700"':'')+' style="font-variant-numeric:tabular-nums">'+(i+1)+' · '+z.rIn+'\u2013'+z.rOut+' · '+(z.bare?'bare':f2(z.local)+' s')+(z.governing?' \u2190 governs':'')+'</text>');
    });
    if(res.zones.length>Math.floor((H-100)/22)) parts.push('<text x="'+LX+'" y="'+(H-54)+'" font-size="10" fill="var(--muted)">\u2026 see the table</text>');
    // colour legend
    var stops=VIRIDIS.map(function(c,i){return '<stop offset="'+(100*i/(VIRIDIS.length-1)).toFixed(1)+'%" stop-color="'+c+'"/>';}).join('');
    parts.push('<defs><linearGradient id="ltGrad" x1="0" x2="1" y1="0" y2="0">'+stops+'</linearGradient></defs>');
    parts.push('<rect x="'+LX+'" y="'+(H-30)+'" width="'+(W-LX-10)+'" height="8" rx="2" fill="url(#ltGrad)"/>');
    parts.push('<text x="'+LX+'" y="'+(H-34)+'" font-size="10" fill="var(--muted)">0 s</text><text x="'+(W-10)+'" y="'+(H-34)+'" text-anchor="end" font-size="10" fill="var(--muted)">'+f2(maxT)+' s</text>');
    parts.push('<text x="'+LX+'" y="'+(H-10)+'" font-size="10" fill="var(--muted)">pulse scale; blank = bare</text>');
    return '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 '+W+' '+H+'" width="'+W+'" height="'+H+'" role="img" aria-label="Wafer map: concentric zones coloured by the pulse the supply rule requires locally, with a zone list and colour scale">'+parts.join('')+'</svg>';
  }

  function init(root){
    root=root||document.getElementById('layout-tool'); if(!root||root.dataset.ltReady) return; root.dataset.ltReady='1';
    var q=function(k){return root.querySelector('[data-lt="'+k+'"]');};
    var tbody=q('zones'), chIn={x:q('x'), Q:q('Q'), p:q('p'), T:q('T'), N0:q('N0')}, csv=q('csv'), csvMsg=q('csvmsg');
    var zones=presets.dense.zones.map(function(z){return Object.assign({},z);});
    var last=null;

    function zoneRows(){ // read the table back into zones
      return Array.prototype.map.call(tbody.querySelectorAll('tr'),function(tr){ var g=function(k){return +tr.querySelector('[data-k="'+k+'"]').value;};
        return {rIn:g('rIn'), rOut:g('rOut'), d:g('d'), AR:g('AR'), density:g('density')}; });
    }
    function cell(k,v,step,min){ return '<td><input data-k="'+k+'" type="number" step="'+step+'" min="'+min+'" value="'+esc(v)+'" aria-label="'+k+'"></td>'; }
    function drawTable(){
      tbody.innerHTML=zones.map(function(z,i){ return '<tr><td class="num lt-idx">'+(i+1)+'</td>'+cell('rIn',z.rIn,'5','0')+cell('rOut',z.rOut,'5','0')+cell('d',z.d,'10','0')+cell('AR',z.AR,'5','0')+
        cell('density',(+z.density||0)?(+z.density).toExponential(2).replace('e+','e'):'0','any','0')+'<td><button type="button" class="lt-x" data-del="'+i+'" aria-label="remove zone '+(i+1)+'" title="remove zone">×</button></td></tr>'; }).join('');
    }
    function chamber(){ return {x:+chIn.x.value, Q:+chIn.Q.value, p:+chIn.p.value, T:+chIn.T.value, N0:+chIn.N0.value*1e18}; }
    function render(){
      zones=zoneRows(); var r=last=compute(zones,chamber()), R=r.recipe, g=r.zones[r.governing];
      q('map').innerHTML=svgMap(r);
      q('recipe').textContent=R.local==null?'—':f2(R.local)+' s';
      q('recipeSub').textContent=g?'local rule of zone '+(r.governing+1)+' ('+g.rIn+'–'+g.rOut+' mm, wall ratio '+f1(g.wall)+'×), conservative':'no patterned zone';
      q('expected').textContent=R.lo==null?'—':rangeTxt(R.lo,R.hi)+' s';
      q('expectedSub').textContent=r.uniform?'expected: 0.98× the rule (uniform wafer, one reference run)':'expected: 0.82–0.87× the local rule (lateral feed; two runs)';
      q('avg').textContent=f2(r.average.rule)+' s'; q('avgSub').textContent='wafer-average rule, wall ratio '+f1(r.average.wall)+'× · under-dose';
      q('under').textContent=r.underDose==null?'—':f2(r.underDose)+'×';
      q('underSub').textContent=r.underDose==null?'under-dose factor':r.underDose<1.05?'average ≈ local: the layout is uniform':'the average rule stops at '+(100/r.underDose).toFixed(0)+' % of what the governing zone needs';
      var sp=r.arSpread, spTxt='';
      if(sp&&sp.ratio>1.001) spTxt='Depth tiers: AR '+sp.max+' to '+sp.min+' ('+((1-sp.min/sp.max)*100).toFixed(0)+' % shallower) gives a '+sp.timeSpreadPct.toFixed(0)+' % spread in the local rule; the chamber run with a 40 % ramp measured 39 %. The shallow tier is over-dosed by about its depth ratio.';
      q('spread').textContent=spTxt;
      var warn=r.warnings.slice(); q('warn').textContent=warn.join(' · ');
      q('results').innerHTML=r.zones.map(function(z,i){
        var exp=z.bare?'<span class="lt-muted">flat '+f2(z.expectedLo)+'–'+f2(z.expectedHi)+'<sup>b</sup></span>':z.governing?'<b>'+rangeTxt(z.expectedLo,z.expectedHi)+'</b>':'<span class="lt-muted">before recipe ends<sup>c</sup></span>';
        return '<tr'+(z.governing?' class="hl"':'')+'><td>'+(i+1)+(z.governing?' <span class="tag">governs</span>':'')+'</td><td class="num">'+z.rIn+'–'+z.rOut+'</td><td>'+(z.bare?'bare':'d '+z.d+' nm · AR '+z.AR+' · '+fden(z.density)+' /m²')+'</td><td class="num">'+f1(z.wall)+'</td><td class="num">'+(z.bare?'<span class="lt-muted">'+f2(z.local)+'<sup>a</sup></span>':f2(z.local))+'</td><td class="num">'+exp+'</td></tr>'; }).join('');
      root.dispatchEvent(new CustomEvent('lt-render',{detail:r}));
    }
    function load(zs){ zones=zs.map(function(z){return Object.assign({},z);}); drawTable(); render(); }

    tbody.addEventListener('input',render);
    tbody.addEventListener('click',function(e){ var b=e.target.closest('[data-del]'); if(!b) return; zones=zoneRows(); zones.splice(+b.dataset.del,1); drawTable(); render(); });
    q('add').addEventListener('click',function(){ zones=zoneRows(); var lastZ=zones[zones.length-1]; var rIn=lastZ?lastZ.rOut:0; zones.push({rIn:rIn, rOut:Math.min(rIn+25,WAFER_R_MM), d:100, AR:100, density:1.78e13}); drawTable(); render(); });
    Object.keys(chIn).forEach(function(k){ chIn[k].addEventListener('input',render); });
    root.querySelectorAll('[data-preset]').forEach(function(b){ b.addEventListener('click',function(){ var p=presets[b.dataset.preset]; if(p){ load(p.zones); csvMsg.textContent=''; } }); });
    q('load').addEventListener('click',function(){ var p=parseCsv(csv.value); if(p.zones.length){ load(p.zones); } csvMsg.textContent=(p.zones.length?p.zones.length+' zones loaded':'no zones read')+(p.errors.length?' · '+p.errors.join('; '):''); });
    q('export').addEventListener('click',function(){ csv.value=toCsv(zoneRows()); csvMsg.textContent='current zones written as CSV'; });
    document.addEventListener('themechange',function(){ if(root.dataset.ltReady) render(); });
    csv.value=toCsv(zones); drawTable(); render();
    return {render:render, load:load, get result(){return last;}};
  }

  window.LayoutTool={compute:compute, parseCsv:parseCsv, toCsv:toCsv, presets:presets, chamberDefault:chamberDefault, init:init, viridis:viridis};
  if(typeof document!=='undefined'){
    var boot=function(){ var el=document.getElementById('layout-tool'); if(el) init(el); else if(document.readyState==='loading') document.addEventListener('DOMContentLoaded',function(){ init(); },{once:true}); };
    if(window.SITE&&window.SITE.ready) boot(); // site.js (defer, earlier in the page) has already dispatched siteready by the time a later defer script runs
    else { document.addEventListener('siteready',boot,{once:true}); if(document.readyState==='loading') document.addEventListener('DOMContentLoaded',function(){ setTimeout(boot,0); },{once:true}); else setTimeout(boot,0); }
  }
})();
