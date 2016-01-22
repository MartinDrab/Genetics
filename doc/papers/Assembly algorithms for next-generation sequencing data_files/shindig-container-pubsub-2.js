var gadgets=gadgets||{};
var shindig=shindig||{};
var osapi=osapi||{};;
gadgets.config=function(){var A={};
var B;
return{register:function(E,D,C){var F=A[E];
if(!F){F=[];
A[E]=F
}F.push({validators:D||{},callback:C})
},get:function(C){if(C){return B[C]||{}
}return B
},init:function(E,L){B=E;
for(var C in A){if(A.hasOwnProperty(C)){var D=A[C],I=E[C];
for(var H=0,G=D.length;
H<G;
++H){var J=D[H];
if(I&&!L){var F=J.validators;
for(var K in F){if(F.hasOwnProperty(K)){if(!F[K](I[K])){throw new Error('Invalid config value "'+I[K]+'" for parameter "'+K+'" in component "'+C+'"')
}}}}if(J.callback){J.callback(E)
}}}}},EnumValidator:function(F){var E=[];
if(arguments.length>1){for(var D=0,C;
(C=arguments[D]);
++D){E.push(C)
}}else{E=F
}return function(H){for(var G=0,I;
(I=E[G]);
++G){if(H===E[G]){return true
}}return false
}
},RegExValidator:function(C){return function(D){return C.test(D)
}
},ExistsValidator:function(C){return typeof C!=="undefined"
},NonEmptyStringValidator:function(C){return typeof C==="string"&&C.length>0
},BooleanValidator:function(C){return typeof C==="boolean"
},LikeValidator:function(C){return function(E){for(var F in C){if(C.hasOwnProperty(F)){var D=C[F];
if(!D(E[F])){return false
}}}return true
}
}}
}();;
gadgets.config.isGadget=false;
gadgets.config.isContainer=true;;
gadgets.util=function(){function G(K){var L;
var I=K.indexOf("?");
var J=K.indexOf("#");
if(J===-1){L=K.substr(I+1)
}else{L=[K.substr(I+1,J-I-1),"&",K.substr(J+1)].join("")
}return L.split("&")
}var E=null;
var D={};
var C={};
var F=[];
var A={0:false,10:true,13:true,34:true,39:true,60:true,62:true,92:true,8232:true,8233:true};
function B(I,J){return String.fromCharCode(J)
}function H(I){D=I["core.util"]||{}
}if(gadgets.config){gadgets.config.register("core.util",null,H)
}return{getUrlParameters:function(R){var J=typeof R==="undefined";
if(E!==null&&J){return E
}var N={};
var K=G(R||document.location.href);
var P=window.decodeURIComponent?decodeURIComponent:unescape;
for(var M=0,L=K.length;
M<L;
++M){var O=K[M].indexOf("=");
if(O===-1){continue
}var I=K[M].substring(0,O);
var Q=K[M].substring(O+1);
Q=Q.replace(/\+/g," ");
N[I]=P(Q)
}if(J){E=N
}return N
},makeClosure:function(L,N,M){var K=[];
for(var J=2,I=arguments.length;
J<I;
++J){K.push(arguments[J])
}return function(){var O=K.slice();
for(var Q=0,P=arguments.length;
Q<P;
++Q){O.push(arguments[Q])
}return N.apply(L,O)
}
},makeEnum:function(J){var K,I,L={};
for(K=0;
(I=J[K]);
++K){L[I]=I
}return L
},getFeatureParameters:function(I){return typeof D[I]==="undefined"?null:D[I]
},hasFeature:function(I){return typeof D[I]!=="undefined"
},getServices:function(){return C
},registerOnLoadHandler:function(I){F.push(I)
},runOnLoadHandlers:function(){for(var J=0,I=F.length;
J<I;
++J){F[J]()
}},escape:function(I,M){if(!I){return I
}else{if(typeof I==="string"){return gadgets.util.escapeString(I)
}else{if(typeof I==="array"){for(var L=0,J=I.length;
L<J;
++L){I[L]=gadgets.util.escape(I[L])
}}else{if(typeof I==="object"&&M){var K={};
for(var N in I){if(I.hasOwnProperty(N)){K[gadgets.util.escapeString(N)]=gadgets.util.escape(I[N],true)
}}return K
}}}}return I
},escapeString:function(M){if(!M){return M
}var J=[],L,N;
for(var K=0,I=M.length;
K<I;
++K){L=M.charCodeAt(K);
N=A[L];
if(N===true){J.push("&#",L,";")
}else{if(N!==false){J.push(M.charAt(K))
}}}return J.join("")
},unescapeString:function(I){if(!I){return I
}return I.replace(/&#([0-9]+);/g,B)
},attachBrowserEvent:function(K,J,L,I){if(typeof K.addEventListener!="undefined"){K.addEventListener(J,L,I)
}else{if(typeof K.attachEvent!="undefined"){K.attachEvent("on"+J,L)
}else{gadgets.warn("cannot attachBrowserEvent: "+J)
}}},removeBrowserEvent:function(K,J,L,I){if(K.removeEventListener){K.removeEventListener(J,L,I)
}else{if(K.detachEvent){K.detachEvent("on"+J,L)
}else{gadgets.warn("cannot removeBrowserEvent: "+J)
}}}}
}();
gadgets.util.getUrlParameters();;
var tamings___=tamings___||[];
tamings___.push(function(A){caja___.whitelistFuncs([[gadgets.util,"escapeString"],[gadgets.util,"getFeatureParameters"],[gadgets.util,"getUrlParameters"],[gadgets.util,"hasFeature"],[gadgets.util,"registerOnLoadHandler"],[gadgets.util,"unescapeString"]])
});;
gadgets.log=(function(){var E=1;
var A=2;
var F=3;
var C=4;
var D=function(I){B(E,I)
};
gadgets.warn=function(I){B(A,I)
};
gadgets.error=function(I){B(F,I)
};
gadgets.setLogLevel=function(I){H=I
};
function B(J,I){if(J<H||!G){return 
}if(J===A&&G.warn){G.warn(I)
}else{if(J===F&&G.error){G.error(I)
}else{if(G.log){G.log(I)
}}}}D.INFO=E;
D.WARNING=A;
D.NONE=C;
var H=E;
var G=window.console?window.console:window.opera?window.opera.postError:undefined;
return D
})();;
var tamings___=tamings___||[];
tamings___.push(function(A){___.grantRead(gadgets.log,"INFO");
___.grantRead(gadgets.log,"WARNING");
___.grantRead(gadgets.log,"ERROR");
___.grantRead(gadgets.log,"NONE");
caja___.whitelistFuncs([[gadgets,"log"],[gadgets,"warn"],[gadgets,"error"],[gadgets,"setLogLevel"]])
});;
if(window.JSON&&window.JSON.parse&&window.JSON.stringify){gadgets.json=(function(){var A=/___$/;
return{parse:function(C){try{return window.JSON.parse(C)
}catch(B){return false
}},stringify:function(C){try{return window.JSON.stringify(C,function(E,D){return !A.test(E)?D:null
})
}catch(B){return null
}}}
})()
}else{gadgets.json=function(){function f(n){return n<10?"0"+n:n
}Date.prototype.toJSON=function(){return[this.getUTCFullYear(),"-",f(this.getUTCMonth()+1),"-",f(this.getUTCDate()),"T",f(this.getUTCHours()),":",f(this.getUTCMinutes()),":",f(this.getUTCSeconds()),"Z"].join("")
};
var m={"\b":"\\b","\t":"\\t","\n":"\\n","\f":"\\f","\r":"\\r",'"':'\\"',"\\":"\\\\"};
function stringify(value){var a,i,k,l,r=/["\\\x00-\x1f\x7f-\x9f]/g,v;
switch(typeof value){case"string":return r.test(value)?'"'+value.replace(r,function(a){var c=m[a];
if(c){return c
}c=a.charCodeAt();
return"\\u00"+Math.floor(c/16).toString(16)+(c%16).toString(16)
})+'"':'"'+value+'"';
case"number":return isFinite(value)?String(value):"null";
case"boolean":case"null":return String(value);
case"object":if(!value){return"null"
}a=[];
if(typeof value.length==="number"&&!value.propertyIsEnumerable("length")){l=value.length;
for(i=0;
i<l;
i+=1){a.push(stringify(value[i])||"null")
}return"["+a.join(",")+"]"
}for(k in value){if(k.match("___$")){continue
}if(value.hasOwnProperty(k)){if(typeof k==="string"){v=stringify(value[k]);
if(v){a.push(stringify(k)+":"+v)
}}}}return"{"+a.join(",")+"}"
}return"undefined"
}return{stringify:stringify,parse:function(text){if(/^[\],:{}\s]*$/.test(text.replace(/\\["\\\/b-u]/g,"@").replace(/"[^"\\\n\r]*"|true|false|null|-?\d+(?:\.\d*)?(?:[eE][+\-]?\d+)?/g,"]").replace(/(?:^|:|,)(?:\s*\[)+/g,""))){return eval("("+text+")")
}return false
}}
}()
}gadgets.json.flatten=function(C){var D={};
if(C===null||C===undefined){return D
}for(var A in C){if(C.hasOwnProperty(A)){var B=C[A];
if(null===B||undefined===B){continue
}D[A]=(typeof B==="string")?B:gadgets.json.stringify(B)
}}return D
};;
var tamings___=tamings___||[];
tamings___.push(function(A){___.tamesTo(gadgets.json.stringify,safeJSON.stringify);
___.tamesTo(gadgets.json.parse,safeJSON.parse)
});;
shindig.Auth=function(){var authToken=null;
var trusted=null;
function addParamsToToken(urlParams){var args=authToken.split("&");
for(var i=0;
i<args.length;
i++){var nameAndValue=args[i].split("=");
if(nameAndValue.length===2){var name=nameAndValue[0];
var value=nameAndValue[1];
if(value==="$"){value=encodeURIComponent(urlParams[name]);
args[i]=name+"="+value
}}}authToken=args.join("&")
}function init(configuration){var urlParams=gadgets.util.getUrlParameters();
var config=configuration["shindig.auth"]||{};
if(config.authToken){authToken=config.authToken
}else{if(urlParams.st){authToken=urlParams.st
}}if(authToken!==null){addParamsToToken(urlParams)
}if(config.trustedJson){trusted=eval("("+config.trustedJson+")")
}}gadgets.config.register("shindig.auth",null,init);
return{getSecurityToken:function(){return authToken
},updateSecurityToken:function(newToken){authToken=newToken
},getTrustedData:function(){return trusted
}}
};;
shindig.auth=new shindig.Auth();;
gadgets.rpctx=gadgets.rpctx||{};
if(!gadgets.rpctx.wpm){gadgets.rpctx.wpm=function(){var F,D;
var C;
var E=false;
var A=false;
var G=false;
function B(){var I=false;
function J(K){if(K.data=="postmessage.test"){I=true;
if(typeof K.origin==="undefined"){A=true
}}}gadgets.util.attachBrowserEvent(window,"message",J,false);
window.postMessage("postmessage.test","*");
if(I){E=true
}gadgets.util.removeBrowserEvent(window,"message",J,false)
}function H(K){var L=gadgets.json.parse(K.data);
if(G){if(!L||!L.f){return 
}var J=gadgets.rpc.getRelayUrl(L.f)||gadgets.util.getUrlParameters()["parent"];
var I=gadgets.rpc.getOrigin(J);
if(!A?K.origin!==I:K.domain!==/^.+:\/\/([^:]+).*/.exec(I)[1]){return 
}}F(L)
}return{getCode:function(){return"wpm"
},isParentVerifiable:function(){return true
},init:function(I,J){F=I;
D=J;
B();
if(!E){C=function(L,M,K){L.postMessage(M,K)
}
}else{C=function(L,M,K){window.setTimeout(function(){L.postMessage(M,K)
},0)
}
}gadgets.util.attachBrowserEvent(window,"message",H,false);
D("..",true);
return true
},setup:function(K,J,I){G=I;
if(K===".."){if(G){gadgets.rpc._createRelayIframe(J)
}else{gadgets.rpc.call(K,gadgets.rpc.ACK)
}}return true
},call:function(J,N,M){var L=gadgets.rpc._getTargetWin(J);
var K=gadgets.rpc.getRelayUrl(J)||gadgets.util.getUrlParameters()["parent"];
var I=gadgets.rpc.getOrigin(K);
if(I){C(L,gadgets.json.stringify(M),I)
}else{gadgets.error("No relay set (used as window.postMessage targetOrigin), cannot send cross-domain message")
}return true
},relayOnload:function(J,I){D(J,true)
}}
}()
};;
gadgets.rpctx=gadgets.rpctx||{};
if(!gadgets.rpctx.frameElement){gadgets.rpctx.frameElement=function(){var E="__g2c_rpc";
var B="__c2g_rpc";
var D;
var C;
function A(G,K,J){try{if(K!==".."){var F=window.frameElement;
if(typeof F[E]==="function"){if(typeof F[E][B]!=="function"){F[E][B]=function(L){D(gadgets.json.parse(L))
}
}F[E](gadgets.json.stringify(J));
return true
}}else{var I=document.getElementById(G);
if(typeof I[E]==="function"&&typeof I[E][B]==="function"){I[E][B](gadgets.json.stringify(J));
return true
}}}catch(H){}return false
}return{getCode:function(){return"fe"
},isParentVerifiable:function(){return false
},init:function(F,G){D=F;
C=G;
return true
},setup:function(J,F){if(J!==".."){try{var I=document.getElementById(J);
I[E]=function(K){D(gadgets.json.parse(K))
}
}catch(H){return false
}}if(J===".."){C("..",true);
var G=function(){window.setTimeout(function(){gadgets.rpc.call(J,gadgets.rpc.ACK)
},500)
};
gadgets.util.registerOnLoadHandler(G)
}return true
},call:function(F,H,G){return A(F,H,G)
}}
}()
};;
gadgets.rpctx=gadgets.rpctx||{};
if(!gadgets.rpctx.nix){gadgets.rpctx.nix=function(){var E="GRPC____NIXVBS_wrapper";
var F="GRPC____NIXVBS_get_wrapper";
var H="GRPC____NIXVBS_handle_message";
var C="GRPC____NIXVBS_create_channel";
var B=10;
var M=500;
var L={};
var D={};
var K;
var J=0;
function G(){var O=L[".."];
if(O){return 
}if(++J>B){gadgets.warn("Nix transport setup failed, falling back...");
K("..",false);
return 
}if(!O&&window.opener&&"GetAuthToken" in window.opener){O=window.opener;
if(O.GetAuthToken()==gadgets.rpc.getAuthToken("..")){var N=gadgets.rpc.getAuthToken("..");
O.CreateChannel(window[F]("..",N),N);
L[".."]=O;
window.opener=null;
K("..",true);
return 
}}window.setTimeout(function(){G()
},M)
}function I(){var O=window.location.href;
var N=O.indexOf("#");
if(N==-1){return O
}return O.substring(0,N)
}function A(P){var O=(2147483647*Math.random())|0;
var Q=[I(),O];
gadgets.rpc._createRelayIframe(P,Q);
var R=window.location.href.split("#")[1]||"";
function N(){var T=window.location.href.split("#")[1]||"";
if(T!==R){clearInterval(S);
var U=gadgets.util.getUrlParameters(window.location.href);
if(U.childtoken==O){G();
return 
}K("..",false)
}}var S=setInterval(N,100)
}return{getCode:function(){return"nix"
},isParentVerifiable:function(N){if(N){return D[N]
}return false
},init:function(O,P){K=P;
if(typeof window[F]!=="unknown"){window[H]=function(R){window.setTimeout(function(){O(gadgets.json.parse(R))
},0)
};
window[C]=function(R,T,S){if(gadgets.rpc.getAuthToken(R)===S){L[R]=T;
K(R,true)
}};
var N="Class "+E+"\n Private m_Intended\nPrivate m_Auth\nPublic Sub SetIntendedName(name)\n If isEmpty(m_Intended) Then\nm_Intended = name\nEnd If\nEnd Sub\nPublic Sub SetAuth(auth)\n If isEmpty(m_Auth) Then\nm_Auth = auth\nEnd If\nEnd Sub\nPublic Sub SendMessage(data)\n "+H+"(data)\nEnd Sub\nPublic Function GetAuthToken()\n GetAuthToken = m_Auth\nEnd Function\nPublic Sub CreateChannel(channel, auth)\n Call "+C+"(m_Intended, channel, auth)\nEnd Sub\nEnd Class\nFunction "+F+"(name, auth)\nDim wrap\nSet wrap = New "+E+"\nwrap.SetIntendedName name\nwrap.SetAuth auth\nSet "+F+" = wrap\nEnd Function";
try{window.execScript(N,"vbscript")
}catch(Q){return false
}}return true
},setup:function(S,O,N){D[S]=!!N;
if(S===".."){if(N){A(O)
}else{G()
}return true
}try{var Q=document.getElementById(S);
var R=window[F](S,O);
Q.contentWindow.opener=R
}catch(P){return false
}return true
},call:function(N,Q,P){try{if(L[N]){L[N].SendMessage(gadgets.json.stringify(P))
}}catch(O){return false
}return true
},relayOnload:function(Q,O){var P=O[0]+"#childtoken="+O[1];
var N=document.getElementById(Q);
N.src=P
}}
}()
};;
gadgets.rpctx=gadgets.rpctx||{};
if(!gadgets.rpctx.flash){gadgets.rpctx.flash=function(){var Z="___xpcswf";
var Q=null;
var J=false;
var K=null;
var g=null;
var U=null;
var h=100;
var R=50;
var a=[];
var i=null;
var A=0;
var V="_scr";
var F="_pnt";
var I=100;
var P=50;
var M=0;
var E=null;
var Y={};
var c=window.location.protocol+"//"+window.location.host;
var N="___jsl";
var D="_fm";
var H;
function T(){window[N]=window[N]||{};
var k=window[N];
var l=k[D]={};
H=N+"."+D;
return l
}var L=T();
function j(m,k){var l=function(){m.apply({},arguments)
};
L[k]=L[k]||l;
return H+"."+k
}function O(k){return k===".."?gadgets.rpc.RPC_ID:k
}function d(k){return k===".."?"INNER":"OUTER"
}function f(k){if(J){Q=k.rpc["commSwf"]||"//xpc.googleusercontent.com/gadgets/xpc.swf"
}}gadgets.config.register("rpc",null,f);
function W(){if(U===null&&document.body&&Q){var m=Q+"?cb="+Math.random()+"&origin="+c+"&jsl=1";
var l=document.createElement("div");
l.style.height="1px";
l.style.width="1px";
var k='<object height="1" width="1" id="'+Z+'" type="application/x-shockwave-flash"><param name="allowScriptAccess" value="always"></param><param name="movie" value="'+m+'"></param><embed type="application/x-shockwave-flash" allowScriptAccess="always" src="'+m+'" height="1" width="1"></embed></object>';
document.body.appendChild(l);
l.innerHTML=k;
U=l.firstChild
}++A;
if(i!==null&&(U!==null||A>=R)){window.clearTimeout(i)
}else{i=window.setTimeout(W,h)
}}function b(){if(Y[".."]){return 
}X("..");
M++;
if(M>=P&&E!==null){window.clearTimeout(E);
E=null
}else{E=window.setTimeout(b,I)
}}function e(){if(U!==null&&U.setup){while(a.length>0){var l=a.shift();
var k=l.targetId;
U.setup(l.token,O(k),d(k))
}}}function G(){e();
if(i!==null){window.clearTimeout(i)
}i=null
}j(G,"ready");
function B(){if(!Y[".."]&&E===null){E=window.setTimeout(b,I)
}}j(B,"setupDone");
function C(m,q,o){var l=gadgets.rpc.getTargetOrigin(m);
var p=gadgets.rpc.getAuthToken(m);
var k="sendMessage_"+O(m)+"_"+p+"_"+d(m);
var n=U[k];
n.call(U,gadgets.json.stringify(o),l);
return true
}function S(m,o,n){var k=gadgets.json.parse(m);
var l=k[V];
if(l){g(l,true);
Y[l]=true;
if(l!==".."){X(l,true)
}return 
}window.setTimeout(function(){K(k,o)
},0)
}j(S,"receiveMessage");
function X(n,m){var k=gadgets.rpc.RPC_ID;
var l={};
l[V]=m?"..":k;
l[F]=k;
C(n,k,l)
}return{getCode:function(){return"flash"
},isParentVerifiable:function(){return true
},init:function(l,k){K=l;
g=k;
J=true;
return true
},setup:function(l,k){a.push({token:k,targetId:l});
if(U===null&&i===null){i=window.setTimeout(W,h)
}return true
},call:C,_receiveMessage:S,_ready:G,_setupDone:B}
}()
};;
gadgets.rpctx=gadgets.rpctx||{};
if(!gadgets.rpctx.rmr){gadgets.rpctx.rmr=function(){var G=500;
var E=10;
var H={};
var B;
var I;
function K(P,N,O,M){var Q=function(){document.body.appendChild(P);
P.src="about:blank";
if(M){P.onload=function(){L(M)
}
}P.src=N+"#"+O
};
if(document.body){Q()
}else{gadgets.util.registerOnLoadHandler(function(){Q()
})
}}function C(O){if(typeof H[O]==="object"){return 
}var P=document.createElement("iframe");
var M=P.style;
M.position="absolute";
M.top="0px";
M.border="0";
M.opacity="0";
M.width="10px";
M.height="1px";
P.id="rmrtransport-"+O;
P.name=P.id;
var N=gadgets.rpc.getRelayUrl(O);
if(!N){N=gadgets.rpc.getOrigin(gadgets.util.getUrlParameters()["parent"])+"/robots.txt"
}H[O]={frame:P,receiveWindow:null,relayUri:N,searchCounter:0,width:10,waiting:true,queue:[],sendId:0,recvId:0};
if(O!==".."){K(P,N,A(O))
}D(O)
}function D(O){var Q=null;
H[O].searchCounter++;
try{var N=gadgets.rpc._getTargetWin(O);
if(O===".."){Q=N.frames["rmrtransport-"+gadgets.rpc.RPC_ID]
}else{Q=N.frames["rmrtransport-.."]
}}catch(P){}var M=false;
if(Q){M=F(O,Q)
}if(!M){if(H[O].searchCounter>E){return 
}window.setTimeout(function(){D(O)
},G)
}}function J(N,P,T,S){var O=null;
if(T!==".."){O=H[".."]
}else{O=H[N]
}if(O){if(P!==gadgets.rpc.ACK){O.queue.push(S)
}if(O.waiting||(O.queue.length===0&&!(P===gadgets.rpc.ACK&&S&&S.ackAlone===true))){return true
}if(O.queue.length>0){O.waiting=true
}var M=O.relayUri+"#"+A(N);
try{O.frame.contentWindow.location=M;
var Q=O.width==10?20:10;
O.frame.style.width=Q+"px";
O.width=Q
}catch(R){return false
}}return true
}function A(N){var O=H[N];
var M={id:O.sendId};
if(O){M.d=Array.prototype.slice.call(O.queue,0);
M.d.push({s:gadgets.rpc.ACK,id:O.recvId})
}return gadgets.json.stringify(M)
}function L(X){var U=H[X];
var Q=U.receiveWindow.location.hash.substring(1);
var Y=gadgets.json.parse(decodeURIComponent(Q))||{};
var N=Y.d||[];
var O=false;
var T=false;
var V=0;
var M=(U.recvId-Y.id);
for(var P=0;
P<N.length;
++P){var S=N[P];
if(S.s===gadgets.rpc.ACK){I(X,true);
if(U.waiting){T=true
}U.waiting=false;
var R=Math.max(0,S.id-U.sendId);
U.queue.splice(0,R);
U.sendId=Math.max(U.sendId,S.id||0);
continue
}O=true;
if(++V<=M){continue
}++U.recvId;
B(S)
}if(O||(T&&U.queue.length>0)){var W=(X==="..")?gadgets.rpc.RPC_ID:"..";
J(X,gadgets.rpc.ACK,W,{ackAlone:O})
}}function F(P,S){var O=H[P];
try{var N=false;
N="document" in S;
if(!N){return false
}N=typeof S.document=="object";
if(!N){return false
}var R=S.location.href;
if(R==="about:blank"){return false
}}catch(M){return false
}O.receiveWindow=S;
function Q(){L(P)
}if(typeof S.attachEvent==="undefined"){S.onresize=Q
}else{S.attachEvent("onresize",Q)
}if(P===".."){K(O.frame,O.relayUri,A(P),P)
}else{L(P)
}return true
}return{getCode:function(){return"rmr"
},isParentVerifiable:function(){return true
},init:function(M,N){B=M;
I=N;
return true
},setup:function(O,M){try{C(O)
}catch(N){gadgets.warn("Caught exception setting up RMR: "+N);
return false
}return true
},call:function(M,O,N){return J(M,N.s,O,N)
}}
}()
};;
gadgets.rpctx=gadgets.rpctx||{};
if(!gadgets.rpctx.ifpc){gadgets.rpctx.ifpc=function(){var H=[];
var E=0;
var D;
var A=2000;
var G={};
function C(K){var I=[];
for(var L=0,J=K.length;
L<J;
++L){I.push(encodeURIComponent(gadgets.json.stringify(K[L])))
}return I.join("&")
}function B(L){var J;
for(var I=H.length-1;
I>=0;
--I){var M=H[I];
try{if(M&&(M.recyclable||M.readyState==="complete")){M.parentNode.removeChild(M);
if(window.ActiveXObject){H[I]=M=null;
H.splice(I,1)
}else{M.recyclable=false;
J=M;
break
}}}catch(K){}}if(!J){J=document.createElement("iframe");
J.style.border=J.style.width=J.style.height="0px";
J.style.visibility="hidden";
J.style.position="absolute";
J.onload=function(){this.recyclable=true
};
H.push(J)
}J.src=L;
window.setTimeout(function(){document.body.appendChild(J)
},0)
}function F(I,K){for(var J=K-1;
J>=0;
--J){if(typeof I[J]==="undefined"){return false
}}return true
}return{getCode:function(){return"ifpc"
},isParentVerifiable:function(){return true
},init:function(I,J){D=J;
D("..",true);
return true
},setup:function(J,I){D(J,true);
return true
},call:function(S,R,Q){var L=gadgets.rpc.getRelayUrl(S);
++E;
if(!L){gadgets.warn("No relay file assigned for IFPC");
return false
}var I=null,J=[];
if(Q.l){var O=Q.a;
I=[L,"#",C([R,E,1,0,C([R,Q.s,"","",R].concat(O))])].join("");
J.push(I)
}else{I=[L,"#",S,"&",R,"@",E,"&"].join("");
var T=encodeURIComponent(gadgets.json.stringify(Q)),N=A-I.length,P=Math.ceil(T.length/N),M=0,K;
while(T.length>0){K=T.substring(0,N);
T=T.substring(N);
J.push([I,P,"&",M,"&",K].join(""));
M+=1
}}do{B(J.shift())
}while(J.length>0);
return true
},_receiveMessage:function(I,N){var O=I[1],M=parseInt(I[2],10),K=parseInt(I[3],10),L=I[I.length-1],J=M===1;
if(M>1){if(!G[O]){G[O]=[]
}G[O][K]=L;
if(F(G[O],M)){L=G[O].join("");
delete G[O];
J=true
}}if(J){N(gadgets.json.parse(decodeURIComponent(L)))
}}}
}()
};;
if(!window.gadgets["rpc"]){gadgets.rpc=function(){var l="__cb";
var t="";
var u="__ack";
var F=500;
var f=10;
var B="|";
var T="callback";
var G="origin";
var R="referer";
var Q={};
var x={};
var c={};
var a={};
var Y=0;
var L={};
var M={};
var r={};
var D={};
var N={};
var d={};
var E=null;
var P=null;
var Z=(window.top!==window.self);
var U=window.name;
var i=function(){};
var p=0;
var z=1;
var A=2;
var W=window.console;
var w=W&&W.log&&function(AF){W.log(AF)
}||function(){};
var s=(function(){function AF(AG){return function(){w(AG+": call ignored")
}
}return{getCode:function(){return"noop"
},isParentVerifiable:function(){return true
},init:AF("init"),setup:AF("setup"),call:AF("call")}
})();
if(gadgets.util){D=gadgets.util.getUrlParameters()
}function j(){if(D.rpctx=="flash"){return gadgets.rpctx.flash
}if(D.rpctx=="rmr"){return gadgets.rpctx.rmr
}return typeof window.postMessage==="function"?gadgets.rpctx.wpm:typeof window.postMessage==="object"?gadgets.rpctx.wpm:window.ActiveXObject?(gadgets.rpctx.ifpc?gadgets.rpctx.ifpc:gadgets.rpctx.nix):navigator.userAgent.indexOf("WebKit")>0?gadgets.rpctx.rmr:navigator.product==="Gecko"?gadgets.rpctx.frameElement:gadgets.rpctx.ifpc
}function K(AK,AI){if(N[AK]){return 
}var AG=g;
if(!AI){AG=s
}N[AK]=AG;
var AF=d[AK]||[];
for(var AH=0;
AH<AF.length;
++AH){var AJ=AF[AH];
AJ.t=e(AK);
AG.call(AK,AJ.f,AJ)
}d[AK]=[]
}var h=false,v=false;
function n(){if(v){return 
}function AF(){h=true
}if(typeof window.addEventListener!="undefined"){window.addEventListener("unload",AF,false)
}else{if(typeof window.attachEvent!="undefined"){window.attachEvent("onunload",AF)
}}v=true
}function J(AF,AJ,AG,AI,AH){if(!a[AJ]||a[AJ]!==AG){gadgets.error("Invalid auth token. "+a[AJ]+" vs "+AG);
i(AJ,A)
}AH.onunload=function(){if(M[AJ]&&!h){i(AJ,z);
gadgets.rpc.removeReceiver(AJ)
}};
n();
AI=gadgets.json.parse(decodeURIComponent(AI))
}function AA(AJ,AG){if(AJ&&typeof AJ.s==="string"&&typeof AJ.f==="string"&&AJ.a instanceof Array){if(a[AJ.f]){if(a[AJ.f]!==AJ.t){gadgets.error("Invalid auth token. "+a[AJ.f]+" vs "+AJ.t);
i(AJ.f,A)
}}if(AJ.s===u){window.setTimeout(function(){K(AJ.f,true)
},0);
return 
}if(AJ.c){AJ[T]=function(AK){gadgets.rpc.call(AJ.f,l,null,AJ.c,AK)
}
}if(AG){var AH=S(AG);
AJ[G]=AG;
var AI=AJ.r;
if(!AI||S(AI)!=AH){AI=AG
}AJ[R]=AI
}var AF=(Q[AJ.s]||Q[t]).apply(AJ,AJ.a);
if(AJ.c&&typeof AF!=="undefined"){gadgets.rpc.call(AJ.f,l,null,AJ.c,AF)
}}}function S(AH){if(!AH){return""
}AH=AH.toLowerCase();
if(AH.indexOf("//")==0){AH=window.location.protocol+AH
}if(AH.indexOf("://")==-1){AH=window.location.protocol+"//"+AH
}var AI=AH.substring(AH.indexOf("://")+3);
var AF=AI.indexOf("/");
if(AF!=-1){AI=AI.substring(0,AF)
}var AK=AH.substring(0,AH.indexOf("://"));
var AJ="";
var AL=AI.indexOf(":");
if(AL!=-1){var AG=AI.substring(AL+1);
AI=AI.substring(0,AL);
if((AK==="http"&&AG!=="80")||(AK==="https"&&AG!=="443")){AJ=":"+AG
}}return AK+"://"+AI+AJ
}function b(AG,AF){return"/"+AG+(AF?B+AF:"")
}function X(AI){if(AI.charAt(0)=="/"){var AG=AI.indexOf(B);
var AH=AG>0?AI.substring(1,AG):AI.substring(1);
var AF=AG>0?AI.substring(AG+1):null;
return{id:AH,origin:AF}
}else{return null
}}function AE(AH){if(typeof AH==="undefined"||AH===".."){return window.parent
}var AG=X(AH);
if(AG){return window.top.frames[AG.id]
}AH=String(AH);
var AF=window.frames[AH];
if(AF){return AF
}AF=document.getElementById(AH);
if(AF&&AF.contentWindow){return AF.contentWindow
}return null
}function k(AI){var AH=null;
var AF=o(AI);
if(AF){AH=AF
}else{var AG=X(AI);
if(AG){AH=AG.origin
}else{if(AI==".."){AH=D.parent
}else{AH=document.getElementById(AI).src
}}}return S(AH)
}var g=j();
Q[t]=function(){w("Unknown RPC service: "+this.s)
};
Q[l]=function(AG,AF){var AH=L[AG];
if(AH){delete L[AG];
AH.call(this,AF)
}};
function y(AH,AF){if(M[AH]===true){return 
}if(typeof M[AH]==="undefined"){M[AH]=0
}var AG=AE(AH);
if(AH===".."||AG!=null){if(g.setup(AH,AF)===true){M[AH]=true;
return 
}}if(M[AH]!==true&&M[AH]++<f){window.setTimeout(function(){y(AH,AF)
},F)
}else{N[AH]=s;
M[AH]=true
}}function m(AG,AJ){if(typeof r[AG]==="undefined"){r[AG]=false;
var AI=o(AG);
if(S(AI)!==S(window.location.href)){return false
}var AH=AE(AG);
try{var AK=AH.gadgets;
r[AG]=AK.rpc.receiveSameDomain
}catch(AF){}}if(typeof r[AG]==="function"){r[AG](AJ);
return true
}return false
}function o(AG){var AF=x[AG];
if(AF&&AF.substring(0,1)==="/"){if(AF.substring(1,2)==="/"){AF=document.location.protocol+AF
}else{AF=document.location.protocol+"//"+document.location.host+AF
}}return AF
}function AD(AG,AF,AH){if(!/http(s)?:\/\/.+/.test(AF)){if(AF.indexOf("//")==0){AF=window.location.protocol+AF
}else{if(AF.charAt(0)=="/"){AF=window.location.protocol+"//"+window.location.host+AF
}else{if(AF.indexOf("://")==-1){AF=window.location.protocol+"//"+AF
}}}}x[AG]=AF;
if(typeof AH!=="undefined"){c[AG]=!!AH
}}function e(AF){return a[AF]
}function C(AF,AG){AG=AG||"";
a[AF]=String(AG);
y(AF,AG)
}function AC(AG){var AF=AG.passReferrer||"";
var AH=AF.split(":",2);
E=AH[0]||"none";
P=AH[1]||"origin"
}function AB(AF){if(q(AF)){g=gadgets.rpctx.ifpc;
g.init(AA,K)
}}function q(AF){return String(AF.useLegacyProtocol)==="true"
}function H(AG,AF){function AH(AK){var AJ=AK?AK.rpc:{};
AC(AJ);
var AI=AJ.parentRelayUrl||"";
AI=S(D.parent||AF)+AI;
AD("..",AI,q(AJ));
AB(AJ);
C("..",AG)
}if(!D.parent&&AF){AH({});
return 
}gadgets.config.register("rpc",null,AH)
}function O(AG,AK,AM){var AJ=null;
if(AG.charAt(0)!="/"){if(!gadgets.util){return 
}AJ=document.getElementById(AG);
if(!AJ){throw new Error("Cannot set up gadgets.rpc receiver with ID: "+AG+", element not found.")
}}var AF=AJ&&AJ.src;
var AH=AK||gadgets.rpc.getOrigin(AF);
AD(AG,AH);
var AL=gadgets.util.getUrlParameters(AF);
var AI=AM||AL.rpctoken;
C(AG,AI)
}function I(AF,AH,AI){if(AF===".."){var AG=AI||D.rpctoken||D.ifpctok||"";
H(AG,AH)
}else{O(AF,AH,AI)
}}function V(AH){if(E==="bidir"||(E==="c2p"&&AH==="..")||(E==="p2c"&&AH!=="..")){var AG=window.location.href;
var AI="?";
if(P==="query"){AI="#"
}else{if(P==="hash"){return AG
}}var AF=AG.lastIndexOf(AI);
AF=AF===-1?AG.length:AF;
return AG.substring(0,AF)
}return null
}return{config:function(AF){if(typeof AF.securityCallback==="function"){i=AF.securityCallback
}},register:function(AG,AF){if(AG===l||AG===u){throw new Error("Cannot overwrite callback/ack service")
}if(AG===t){throw new Error("Cannot overwrite default service: use registerDefault")
}Q[AG]=AF
},unregister:function(AF){if(AF===l||AF===u){throw new Error("Cannot delete callback/ack service")
}if(AF===t){throw new Error("Cannot delete default service: use unregisterDefault")
}delete Q[AF]
},registerDefault:function(AF){Q[t]=AF
},unregisterDefault:function(){delete Q[t]
},forceParentVerifiable:function(){if(!g.isParentVerifiable()){g=gadgets.rpctx.ifpc
}},call:function(AF,AH,AM,AK){AF=AF||"..";
var AL="..";
if(AF===".."){AL=U
}else{if(AF.charAt(0)=="/"){AL=b(U,gadgets.rpc.getOrigin(window.location.href))
}}++Y;
if(AM){L[Y]=AM
}var AJ={s:AH,f:AL,c:AM?Y:0,a:Array.prototype.slice.call(arguments,3),t:a[AF],l:!!c[AF]};
var AG=V(AF);
if(AG){AJ.r=AG
}if(AF!==".."&&X(AF)==null&&!document.getElementById(AF)){return 
}if(m(AF,AJ)){return 
}var AI=N[AF];
if(!AI&&X(AF)!==null){AI=g
}if(!AI){if(!d[AF]){d[AF]=[AJ]
}else{d[AF].push(AJ)
}return 
}if(c[AF]){AI=gadgets.rpctx.ifpc
}if(AI.call(AF,AL,AJ)===false){N[AF]=s;
g.call(AF,AL,AJ)
}},getRelayUrl:o,setRelayUrl:AD,setAuthToken:C,setupReceiver:I,getAuthToken:e,removeReceiver:function(AF){delete x[AF];
delete c[AF];
delete a[AF];
delete M[AF];
delete r[AF];
delete N[AF]
},getRelayChannel:function(){return g.getCode()
},receive:function(AG,AF){if(AG.length>4){g._receiveMessage(AG,AA)
}else{J.apply(null,AG.concat(AF))
}},receiveSameDomain:function(AF){AF.a=Array.prototype.slice.call(AF.a);
window.setTimeout(function(){AA(AF)
},0)
},getOrigin:S,getTargetOrigin:k,init:function(){if(g.init(AA,K)===false){g=s
}if(Z){I("..")
}else{gadgets.config.register("rpc",null,function(AG){var AF=AG.rpc||{};
AC(AF);
AB(AF)
})
}},_getTargetWin:AE,_parseSiblingId:X,ACK:u,RPC_ID:U||"..",SEC_ERROR_LOAD_TIMEOUT:p,SEC_ERROR_FRAME_PHISH:z,SEC_ERROR_FORGED_MSG:A}
}();
gadgets.rpc.init()
};;
gadgets.io=function(){var config={};
var oauthState;
function makeXhr(){var x;
if(typeof shindig!="undefined"&&shindig.xhrwrapper&&shindig.xhrwrapper.createXHR){return shindig.xhrwrapper.createXHR()
}else{if(typeof ActiveXObject!="undefined"){x=new ActiveXObject("Msxml2.XMLHTTP");
if(!x){x=new ActiveXObject("Microsoft.XMLHTTP")
}return x
}else{if(typeof XMLHttpRequest!="undefined"||window.XMLHttpRequest){return new window.XMLHttpRequest()
}else{throw ("no xhr available")
}}}}function hadError(xobj,callback){if(xobj.readyState!==4){return true
}try{if(xobj.status!==200){var error=(""+xobj.status);
if(xobj.responseText){error=error+" "+xobj.responseText
}callback({errors:[error],rc:xobj.status,text:xobj.responseText});
return true
}}catch(e){callback({errors:[e.number+" Error not specified"],rc:e.number,text:e.description});
return true
}return false
}function processNonProxiedResponse(url,callback,params,xobj){if(hadError(xobj,callback)){return 
}var data={body:xobj.responseText};
callback(transformResponseData(params,data))
}var UNPARSEABLE_CRUFT="throw 1; < don't be evil' >";
function processResponse(url,callback,params,xobj){if(hadError(xobj,callback)){return 
}var txt=xobj.responseText;
var offset=txt.indexOf(UNPARSEABLE_CRUFT)+UNPARSEABLE_CRUFT.length;
if(offset<UNPARSEABLE_CRUFT.length){return 
}txt=txt.substr(offset);
var data=eval("("+txt+")");
data=data[url];
if(data.oauthState){oauthState=data.oauthState
}if(data.st){shindig.auth.updateSecurityToken(data.st)
}callback(transformResponseData(params,data))
}function transformResponseData(params,data){var resp={text:data.body,rc:data.rc||200,headers:data.headers,oauthApprovalUrl:data.oauthApprovalUrl,oauthError:data.oauthError,oauthErrorText:data.oauthErrorText,errors:[]};
if(resp.rc<200||resp.rc>=400){resp.errors=[resp.rc+" Error"]
}else{if(resp.text){if(resp.rc>=300&&resp.rc<400){params.CONTENT_TYPE="TEXT"
}switch(params.CONTENT_TYPE){case"JSON":case"FEED":resp.data=gadgets.json.parse(resp.text);
if(!resp.data){resp.errors.push("500 Failed to parse JSON");
resp.rc=500;
resp.data=null
}break;
case"DOM":var dom;
if(typeof ActiveXObject!="undefined"){dom=new ActiveXObject("Microsoft.XMLDOM");
dom.async=false;
dom.validateOnParse=false;
dom.resolveExternals=false;
if(!dom.loadXML(resp.text)){resp.errors.push("500 Failed to parse XML");
resp.rc=500
}else{resp.data=dom
}}else{var parser=new DOMParser();
dom=parser.parseFromString(resp.text,"text/xml");
if("parsererror"===dom.documentElement.nodeName){resp.errors.push("500 Failed to parse XML");
resp.rc=500
}else{resp.data=dom
}}break;
default:resp.data=resp.text;
break
}}}return resp
}function makeXhrRequest(realUrl,proxyUrl,callback,paramData,method,params,processResponseFunction,opt_contentType){var xhr=makeXhr();
if(proxyUrl.indexOf("//")==0){proxyUrl=document.location.protocol+proxyUrl
}xhr.open(method,proxyUrl,true);
if(callback){xhr.onreadystatechange=gadgets.util.makeClosure(null,processResponseFunction,realUrl,callback,params,xhr)
}if(paramData!==null){xhr.setRequestHeader("Content-Type",opt_contentType||"application/x-www-form-urlencoded");
xhr.send(paramData)
}else{xhr.send(null)
}}function respondWithPreload(postData,params,callback){if(gadgets.io.preloaded_&&postData.httpMethod==="GET"){for(var i=0;
i<gadgets.io.preloaded_.length;
i++){var preload=gadgets.io.preloaded_[i];
if(preload&&(preload.id===postData.url)){delete gadgets.io.preloaded_[i];
if(preload.rc!==200){callback({rc:preload.rc,errors:[preload.rc+" Error"]})
}else{if(preload.oauthState){oauthState=preload.oauthState
}var resp={body:preload.body,rc:preload.rc,headers:preload.headers,oauthApprovalUrl:preload.oauthApprovalUrl,oauthError:preload.oauthError,oauthErrorText:preload.oauthErrorText,errors:[]};
callback(transformResponseData(params,resp))
}return true
}}}return false
}function init(configuration){config=configuration["core.io"]||{}
}var requiredConfig={proxyUrl:new gadgets.config.RegExValidator(/.*%(raw)?url%.*/),jsonProxyUrl:gadgets.config.NonEmptyStringValidator};
gadgets.config.register("core.io",requiredConfig,init);
return{makeRequest:function(url,callback,opt_params){var params=opt_params||{};
var httpMethod=params.METHOD||"GET";
var refreshInterval=params.REFRESH_INTERVAL;
var auth,st;
if(params.AUTHORIZATION&&params.AUTHORIZATION!=="NONE"){auth=params.AUTHORIZATION.toLowerCase();
st=shindig.auth.getSecurityToken()
}else{if(httpMethod==="GET"&&refreshInterval===undefined){refreshInterval=3600
}}var signOwner=true;
if(typeof params.OWNER_SIGNED!=="undefined"){signOwner=params.OWNER_SIGNED
}var signViewer=true;
if(typeof params.VIEWER_SIGNED!=="undefined"){signViewer=params.VIEWER_SIGNED
}var headers=params.HEADERS||{};
if(httpMethod==="POST"&&!headers["Content-Type"]){headers["Content-Type"]="application/x-www-form-urlencoded"
}var urlParams=gadgets.util.getUrlParameters();
var paramData={url:url,httpMethod:httpMethod,headers:gadgets.io.encodeValues(headers,false),postData:params.POST_DATA||"",authz:auth||"",st:st||"",contentType:params.CONTENT_TYPE||"TEXT",numEntries:params.NUM_ENTRIES||"3",getSummaries:!!params.GET_SUMMARIES,signOwner:signOwner,signViewer:signViewer,gadget:urlParams.url,container:urlParams.container||urlParams.synd||"default",bypassSpecCache:gadgets.util.getUrlParameters().nocache||"",nocache:gadgets.util.getUrlParameters(url).nocache||"",getFullHeaders:!!params.GET_FULL_HEADERS};
if(auth==="oauth"||auth==="signed"){if(gadgets.io.oauthReceivedCallbackUrl_){paramData.OAUTH_RECEIVED_CALLBACK=gadgets.io.oauthReceivedCallbackUrl_;
gadgets.io.oauthReceivedCallbackUrl_=null
}paramData.oauthState=oauthState||"";
for(var opt in params){if(params.hasOwnProperty(opt)){if(opt.indexOf("OAUTH_")===0){paramData[opt]=params[opt]
}}}}var proxyUrl=config.jsonProxyUrl.replace("%host%",document.location.host);
if(!respondWithPreload(paramData,params,callback,processResponse)){if(httpMethod==="GET"&&refreshInterval>0){var extraparams="?refresh="+refreshInterval+"&"+gadgets.io.encodeValues(paramData);
makeXhrRequest(url,proxyUrl+extraparams,callback,null,"GET",params,processResponse)
}else{makeXhrRequest(url,proxyUrl,callback,gadgets.io.encodeValues(paramData),"POST",params,processResponse)
}}},makeNonProxiedRequest:function(relativeUrl,callback,opt_params,opt_contentType){var params=opt_params||{};
makeXhrRequest(relativeUrl,relativeUrl,callback,params.POST_DATA,params.METHOD,params,processNonProxiedResponse,opt_contentType)
},clearOAuthState:function(){oauthState=undefined
},encodeValues:function(fields,opt_noEscaping){var escape=!opt_noEscaping;
var buf=[];
var first=false;
for(var i in fields){if(fields.hasOwnProperty(i)&&!/___$/.test(i)){if(!first){first=true
}else{buf.push("&")
}buf.push(escape?encodeURIComponent(i):i);
buf.push("=");
buf.push(escape?encodeURIComponent(fields[i]):fields[i])
}}return buf.join("")
},getProxyUrl:function(url,opt_params){var params=opt_params||{};
var refresh=params.REFRESH_INTERVAL;
if(refresh===undefined){refresh="3600"
}var urlParams=gadgets.util.getUrlParameters();
var rewriteMimeParam=params.rewriteMime?"&rewriteMime="+encodeURIComponent(params.rewriteMime):"";
var ret=config.proxyUrl.replace("%url%",encodeURIComponent(url)).replace("%host%",document.location.host).replace("%rawurl%",url).replace("%refresh%",encodeURIComponent(refresh)).replace("%gadget%",encodeURIComponent(urlParams.url)).replace("%container%",encodeURIComponent(urlParams.container||urlParams.synd||"default")).replace("%rewriteMime%",rewriteMimeParam);
if(ret.indexOf("//")==0){ret=window.location.protocol+ret
}return ret
}}
}();
gadgets.io.RequestParameters=gadgets.util.makeEnum(["METHOD","CONTENT_TYPE","POST_DATA","HEADERS","AUTHORIZATION","NUM_ENTRIES","GET_SUMMARIES","GET_FULL_HEADERS","REFRESH_INTERVAL","OAUTH_SERVICE_NAME","OAUTH_USE_TOKEN","OAUTH_TOKEN_NAME","OAUTH_REQUEST_TOKEN","OAUTH_REQUEST_TOKEN_SECRET","OAUTH_RECEIVED_CALLBACK"]);
gadgets.io.MethodType=gadgets.util.makeEnum(["GET","POST","PUT","DELETE","HEAD"]);
gadgets.io.ContentType=gadgets.util.makeEnum(["TEXT","DOM","JSON","FEED"]);
gadgets.io.AuthorizationType=gadgets.util.makeEnum(["NONE","SIGNED","OAUTH"]);;
var tamings___=tamings___||[];
tamings___.push(function(A){caja___.whitelistFuncs([[gadgets.io,"encodeValues"],[gadgets.io,"getProxyUrl"],[gadgets.io,"makeRequest"]])
});;
shindig.uri=(function(){var A=new RegExp("^(?:([^:/?#]+):)?(?://([^/?#]*))?([^?#]*)(?:\\?([^#]*))?(?:#(.*))?");
return function(Y){var R="";
var N="";
var C="";
var H="";
var D=null;
var I="";
var J=null;
var L=window.decodeURIComponent?decodeURIComponent:unescape;
var X=window.encodeURIComponent?encodeURIComponent:escape;
var K=null;
function U(a){if(a.match(A)===null){throw"Malformed URL: "+a
}R=RegExp.$1;
N=RegExp.$2;
C=RegExp.$3;
H=RegExp.$4;
I=RegExp.$5
}function T(f){var e=[];
for(var c=0,a=f.length;
c<a;
++c){var b=f[c][0];
var d=f[c][1];
if(d===undefined){continue
}e.push(X(b)+(d!==null?"="+X(d):""))
}return e.join("&")
}function Q(){if(D){H=T(D);
D=null
}return H
}function Z(){if(J){I=T(J);
J=null
}return I
}function O(a){D=D||F(H);
return S(D,a)
}function W(a){J=J||F(I);
return S(J,a)
}function B(b,a){D=M(D||F(H),b,a);
return K
}function G(b,a){J=M(J||F(I),b,a);
return K
}function V(){return[R,R!==""?":":"",N!==""?"//":"",N].join("")
}function P(){var b=Q();
var a=Z();
return[V(),C,b!==""?"?":"",b,a!==""?"#":"",a].join("")
}function F(h){var g=[];
var f=h.split("&");
for(var c=0,a=f.length;
c<a;
++c){var e=f[c].split("=");
var b=e.shift();
var d=null;
if(e.length>0){d=e.join("").replace(/\+/g," ")
}g.push([b,d!=null?L(d):null])
}return g
}function S(a,d){for(var c=0,b=a.length;
c<b;
++c){if(a[c][0]==d){return a[c][1]
}}return undefined
}function M(e,f,d){var h=f;
if(typeof f==="string"){h={};
h[f]=d
}for(var c in h){var g=false;
for(var b=0,a=e.length;
!g&&b<a;
++b){if(e[b][0]==c){e[b][1]=h[c];
g=true
}}if(!g){e.push([c,h[c]])
}}return e
}function E(b,a){b=b||"";
if(b[0]===a){b=b.substr(a.length)
}return b
}if(typeof Y==="object"&&typeof Y.toString==="function"){U(Y.toString())
}else{if(Y){U(Y)
}}K={getSchema:function(){return R
},getAuthority:function(){return N
},getOrigin:V,getPath:function(){return C
},getQuery:Q,getFragment:Z,getQP:O,getFP:W,setSchema:function(a){R=a;
return K
},setAuthority:function(a){N=a;
return K
},setPath:function(a){C=(a[0]==="/"?"":"/")+a;
return K
},setQuery:function(a){D=null;
H=E(a,"?");
return K
},setFragment:function(a){J=null;
I=E(a,"#");
return K
},setQP:B,setFP:G,setExistingP:function(a,b){if(O(a,b)!==undefined){B(a,b)
}if(W(a,b)!==undefined){G(a,b)
}return K
},toString:P};
return K
}
})();;
(function(){osapi._registerMethod=function(G,F){var A=typeof ___!=="undefined";
if(G=="newBatch"){return 
}var D=G.split(".");
var C=osapi;
for(var B=0;
B<D.length-1;
B++){C[D[B]]=C[D[B]]||{};
C=C[D[B]]
}var E=function(J){var I=osapi.newBatch();
var H={};
H.execute=function(M){var K=A?___.untame(M):M;
var L=A?___.USELESS:this;
I.add(G,this);
I.execute(function(N){if(N.error){K.call(L,N.error)
}else{K.call(L,N[G])
}})
};
if(A){___.markInnocent(H.execute,"execute")
}J=J||{};
J.userId=J.userId||"@viewer";
J.groupId=J.groupId||"@self";
H.method=G;
H.transport=F;
H.rpc=J;
return H
};
if(A&&typeof ___.markInnocent!=="undefined"){___.markInnocent(E,G)
}if(C[D[D.length-1]]){gadgets.warn("Skipping duplicate osapi method definition "+G+" on transport "+F.name)
}else{C[D[D.length-1]]=E
}}
})();;
(function(){var A=function(){var C={};
var B=[];
var F=function(G,H){if(H&&G){B.push({key:G,request:H})
}return C
};
var E=function(H){var G={method:H.request.method,id:H.key};
if(H.request.rpc){G.params=H.request.rpc
}return G
};
var D=function(G){var H={};
var O={};
var J=0;
var K=[];
for(var M=0;
M<B.length;
M++){var I=B[M].request.transport;
if(!O[I.name]){K.push(I);
J++
}O[I.name]=O[I.name]||[];
O[I.name].push(E(B[M]))
}var N=function(S){if(S.error){H.error=S.error
}for(var R=0;
R<B.length;
R++){var Q=B[R].key;
var P=S[Q];
if(P){if(P.error){H[Q]=P
}else{H[Q]=P.data||P.result
}}}J--;
if(J===0){G(H)
}};
for(var L=0;
L<K.length;
L++){K[L].execute(O[K[L].name],N)
}if(J==0){window.setTimeout(function(){G(H)
},0)
}};
C.execute=D;
C.add=F;
return C
};
osapi.newBatch=A
})();;
(function(){function A(H,G){function F(J){if(J.errors[0]){G({error:{code:J.rc,message:J.text}})
}else{var K=J.result||J.data;
if(K.error){G(K)
}else{var I={};
for(var L=0;
L<K.length;
L++){I[K[L].id]=K[L]
}G(I)
}}}var E={POST_DATA:gadgets.json.stringify(H),CONTENT_TYPE:"JSON",METHOD:"POST",AUTHORIZATION:"SIGNED"};
var C=this.name;
var D=shindig.auth.getSecurityToken();
if(D){C+="?st=";
C+=encodeURIComponent(D)
}gadgets.io.makeNonProxiedRequest(C,F,E,"application/json")
}function B(F){var H=F["osapi.services"];
if(H){for(var E in H){if(H.hasOwnProperty(E)){if(E.indexOf("http")==0||E.indexOf("//")==0){var C=E.replace("%host%",document.location.host);
var I={name:C,execute:A};
var D=H[E];
for(var G=0;
G<D.length;
G++){osapi._registerMethod(D[G],I)
}}}}}}if(gadgets.config){gadgets.config.register("osapi.services",null,B)
}})();;
if(gadgets&&gadgets.rpc){(function(){function A(E,D){var C=function(G){if(!G){D({code:500,message:"Container refused the request"})
}else{if(G.error){D(G)
}else{var F={};
for(var H=0;
H<G.length;
H++){F[G[H].id]=G[H]
}D(F)
}}};
gadgets.rpc.call("..","osapi._handleGadgetRpcMethod",C,E)
}function B(C){var F={name:"gadgets.rpc",execute:A};
var K=C["osapi.services"];
if(K){for(var D in K){if(K.hasOwnProperty(D)){if(D==="gadgets.rpc"){var E=K[D];
for(var H=0;
H<E.length;
H++){osapi._registerMethod(E[H],F)
}}}}}if(osapi.container&&osapi.container.listMethods){var G=gadgets.util.runOnLoadHandlers;
var I=2;
var J=function(){I--;
if(I==0){G()
}};
gadgets.util.runOnLoadHandlers=J;
osapi.container.listMethods({}).execute(function(L){if(!L.error){for(var M=0;
M<L.length;
M++){if(L[M]!="container.listMethods"){osapi._registerMethod(L[M],F)
}}}J()
});
window.setTimeout(J,500)
}}if(gadgets.config&&gadgets.config.isGadget){gadgets.config.register("osapi.services",null,B)
}})()
};;
gadgets.util.registerOnLoadHandler(function(){if(osapi&&osapi.people&&osapi.people.get){osapi.people.getViewer=function(A){A=A||{};
A.userId="@viewer";
A.groupId="@self";
return osapi.people.get(A)
};
osapi.people.getViewerFriends=function(A){A=A||{};
A.userId="@viewer";
A.groupId="@friends";
return osapi.people.get(A)
};
osapi.people.getOwner=function(A){A=A||{};
A.userId="@owner";
A.groupId="@self";
return osapi.people.get(A)
};
osapi.people.getOwnerFriends=function(A){A=A||{};
A.userId="@owner";
A.groupId="@friends";
return osapi.people.get(A)
}
}});;
var tamings___=tamings___||[];
tamings___.push(function(A){___.tamesTo(osapi.newBatch,___.markFuncFreeze(function(){var C=osapi.newBatch();
___.markInnocent(C.add,"add");
___.markInnocent(C.execute,"execute");
return ___.tame(C)
}));
A.outers.osapi=___.tame(osapi);
___.grantRead(A.outers,"osapi");
var B=A;
gadgets.util.registerOnLoadHandler(function(){if(osapi&&osapi.people&&osapi.people.get){caja___.whitelistFuncs([[osapi.people,"getViewer"],[osapi.people,"getViewerFriends"],[osapi.people,"getOwner"],[osapi.people,"getOwnerFriends"]]);
B.outers.osapi.people.getViewer=___.tame(osapi.people.getViewer);
B.outers.osapi.people.getViewerFriends=___.tame(osapi.people.getViewerFriends);
B.outers.osapi.people.getOwner=___.tame(osapi.people.getOwner);
B.outers.osapi.people.getOwnerFriends=___.tame(osapi.people.getOwnerFriends)
}})
});;
shindig._uri=shindig.uri;
shindig.uri=(function(){var C=shindig._uri;
shindig._uri=null;
function A(E,D){return E.getOrigin()==D.getOrigin()
}function B(E,G){if(E.getSchema()==""){E.setSchema(G.getSchema())
}if(E.getAuthority()==""){E.setAuthority(G.getAuthority())
}var F=E.getPath();
if(F==""||F.charAt(0)!="/"){var H=G.getPath();
var D=H.lastIndexOf("/");
if(D!=-1){H=H.substring(0,D+1)
}E.setPath(G.getPath()+F)
}}return function(D){var E=C(D);
E.hasSameOrigin=function(F){return A(E,F)
};
E.resolve=function(F){return B(E,F)
};
return E
}
})();;
Function.prototype.inherits=function(A){function B(){}B.prototype=A.prototype;
this.superClass_=A.prototype;
this.prototype=new B();
this.prototype.constructor=this
};;
shindig.cookies={};
shindig.cookies.JsType_={UNDEFINED:"undefined"};
shindig.cookies.isDef=function(A){return typeof A!=shindig.cookies.JsType_.UNDEFINED
};
shindig.cookies.set=function(C,I,H,B,D){if(/;=/g.test(C)){throw new Error('Invalid cookie name "'+C+'"')
}if(/;/g.test(I)){throw new Error('Invalid cookie value "'+I+'"')
}if(!shindig.cookies.isDef(H)){H=-1
}var F=D?";domain="+D:"";
var J=B?";path="+B:"";
var E;
if(H<0){E=""
}else{if(H===0){var G=new Date(1970,1,1);
E=";expires="+G.toUTCString()
}else{var A=new Date((new Date).getTime()+H*1000);
E=";expires="+A.toUTCString()
}}document.cookie=C+"="+I+F+J+E
};
shindig.cookies.get=function(B,G){var F=B+"=";
var D=String(document.cookie);
for(var H=-1;
(H=D.indexOf(F,H+1))>=0;
){var C=H;
while(--C>=0){var E=D.charAt(C);
if(E==";"){C=-1;
break
}}if(C==-1){var A=D.indexOf(";",H);
if(A<0){A=D.length
}return D.substring(H+F.length,A)
}}return G
};
shindig.cookies.remove=function(B,A,C){var D=shindig.cookies.containsKey(B);
shindig.cookies.set(B,"",0,A,C);
return D
};
shindig.cookies.getKeyValues_=function(){var E=String(document.cookie);
var G=E.split(/\s*;\s*/);
var F=[],A=[],C,B;
for(var D=0;
B=G[D];
D++){C=B.indexOf("=");
if(C==-1){F.push("");
A.push(B)
}else{F.push(B.substring(0,C));
A.push(B.substring(C+1))
}}return{keys:F,values:A}
};
shindig.cookies.getKeys=function(){return shindig.cookies.getKeyValues_().keys
};
shindig.cookies.getValues=function(){return shindig.cookies.getKeyValues_().values
};
shindig.cookies.isEmpty=function(){return document.cookie===""
};
shindig.cookies.getCount=function(){var A=String(document.cookie);
if(A===""){return 0
}var B=A.split(/\s*;\s*/);
return B.length
};
shindig.cookies.containsKey=function(B){var A={};
return shindig.cookies.get(B,A)!==A
};
shindig.cookies.containsValue=function(C){var A=shindig.cookies.getKeyValues_().values;
for(var B=0;
B<A.length;
B++){if(A[B]==C){return true
}}return false
};
shindig.cookies.clear=function(){var B=shindig.cookies.getKeyValues_().keys;
for(var A=B.length-1;
A>=0;
A--){shindig.cookies.remove(B[A])
}};
shindig.cookies.MAX_COOKIE_LENGTH=3950;;
shindig.errors={};
shindig.errors.SUBCLASS_RESPONSIBILITY="subclass responsibility";
shindig.errors.TO_BE_DONE="to be done";
shindig.callAsyncAndJoin=function(E,A,D){var F=E.length;
var C=[];
for(var B=0;
B<E.length;
B++){var G=function(H){var I=E[H];
if(typeof I==="string"){I=D[I]
}I.call(D,function(J){C[H]=J;
if(--F===0){A(C)
}})
};
G(B)
}};
shindig.Extensible=function(){};
shindig.Extensible.prototype.setDependencies=function(A){for(var B in A){this[B]=A[B]
}};
shindig.Extensible.prototype.getDependencies=function(A){return this[A]
};
shindig.UserPrefStore=function(){};
shindig.UserPrefStore.prototype.getPrefs=function(A){throw Error(shindig.errors.SUBCLASS_RESPONSIBILITY)
};
shindig.UserPrefStore.prototype.savePrefs=function(A){throw Error(shindig.errors.SUBCLASS_RESPONSIBILITY)
};
shindig.DefaultUserPrefStore=function(){shindig.UserPrefStore.call(this)
};
shindig.DefaultUserPrefStore.inherits(shindig.UserPrefStore);
shindig.DefaultUserPrefStore.prototype.getPrefs=function(A){};
shindig.DefaultUserPrefStore.prototype.savePrefs=function(A){};
shindig.GadgetService=function(){};
shindig.GadgetService.prototype.setHeight=function(B,A){throw Error(shindig.errors.SUBCLASS_RESPONSIBILITY)
};
shindig.GadgetService.prototype.setTitle=function(A,B){throw Error(shindig.errors.SUBCLASS_RESPONSIBILITY)
};
shindig.GadgetService.prototype.setUserPref=function(A){throw Error(shindig.errors.SUBCLASS_RESPONSIBILITY)
};
shindig.IfrGadgetService=function(){shindig.GadgetService.call(this);
gadgets.rpc.register("resize_iframe",this.setHeight);
gadgets.rpc.register("set_pref",this.setUserPref);
gadgets.rpc.register("set_title",this.setTitle);
gadgets.rpc.register("requestNavigateTo",this.requestNavigateTo);
gadgets.rpc.register("requestSendMessage",this.requestSendMessage)
};
shindig.IfrGadgetService.inherits(shindig.GadgetService);
shindig.IfrGadgetService.prototype.setHeight=function(A){if(A>shindig.container.maxheight_){A=shindig.container.maxheight_
}var B=document.getElementById(this.f);
if(B){B.style.height=A+"px"
}};
shindig.IfrGadgetService.prototype.setTitle=function(B){var A=document.getElementById(this.f+"_title");
if(A){A.innerHTML=B.replace(/&/g,"&amp;").replace(/</g,"&lt;")
}};
shindig.IfrGadgetService.prototype.setUserPref=function(G,B,D){var F=shindig.container.gadgetService.getGadgetIdFromModuleId(this.f);
var E=shindig.container.getGadget(F);
for(var C=1,A=arguments.length;
C<A;
C+=2){this.userPrefs[arguments[C]].value=arguments[C+1]
}E.saveUserPrefs()
};
shindig.IfrGadgetService.prototype.requestSendMessage=function(A,D,B,C){if(B){window.setTimeout(function(){B(new opensocial.ResponseItem(null,null,opensocial.ResponseItem.Error.NOT_IMPLEMENTED,null))
},0)
}};
shindig.IfrGadgetService.prototype.requestNavigateTo=function(A,D){var E=shindig.container.gadgetService.getGadgetIdFromModuleId(this.f);
var B=shindig.container.gadgetService.getUrlForView(A);
if(D){var C=gadgets.json.stringify(D);
if(C.length>0){B+="&appParams="+encodeURIComponent(C)
}}if(B&&document.location.href.indexOf(B)==-1){document.location.href=B
}};
shindig.IfrGadgetService.prototype.getUrlForView=function(A){if(A==="canvas"){return"/canvas"
}else{if(A==="profile"){return"/profile"
}else{return null
}}};
shindig.IfrGadgetService.prototype.getGadgetIdFromModuleId=function(A){return parseInt(A.match(/_([0-9]+)$/)[1],10)
};
shindig.LayoutManager=function(){};
shindig.LayoutManager.prototype.getGadgetChrome=function(A){throw Error(shindig.errors.SUBCLASS_RESPONSIBILITY)
};
shindig.StaticLayoutManager=function(){shindig.LayoutManager.call(this)
};
shindig.StaticLayoutManager.inherits(shindig.LayoutManager);
shindig.StaticLayoutManager.prototype.setGadgetChromeIds=function(A){this.gadgetChromeIds_=A
};
shindig.StaticLayoutManager.prototype.getGadgetChrome=function(B){var A=this.gadgetChromeIds_[B.id];
return A?document.getElementById(A):null
};
shindig.FloatLeftLayoutManager=function(A){shindig.LayoutManager.call(this);
this.layoutRootId_=A
};
shindig.FloatLeftLayoutManager.inherits(shindig.LayoutManager);
shindig.FloatLeftLayoutManager.prototype.getGadgetChrome=function(C){var B=document.getElementById(this.layoutRootId_);
if(B){var A=document.createElement("div");
A.className="gadgets-gadget-chrome";
A.style.cssFloat="left";
B.appendChild(A);
return A
}else{return null
}};
shindig.Gadget=function(B){this.userPrefs={};
if(B){for(var A in B){if(B.hasOwnProperty(A)){this[A]=B[A]
}}}if(!this.secureToken){this.secureToken="john.doe:john.doe:appid:cont:url:0:default"
}};
shindig.Gadget.prototype.getUserPrefs=function(){return this.userPrefs
};
shindig.Gadget.prototype.saveUserPrefs=function(){shindig.container.userPrefStore.savePrefs(this)
};
shindig.Gadget.prototype.getUserPrefValue=function(B){var A=this.userPrefs[B];
return typeof (A.value)!="undefined"&&A.value!=null?A.value:A["default"]
};
shindig.Gadget.prototype.render=function(A){if(A){var B=this;
this.getContent(function(C){A.innerHTML=C;
B.finishRender(A)
})
}};
shindig.Gadget.prototype.getContent=function(A){shindig.callAsyncAndJoin(["getTitleBarContent","getUserPrefsDialogContent","getMainContent"],function(B){A(B.join(""))
},this)
};
shindig.Gadget.prototype.getTitleBarContent=function(A){throw Error(shindig.errors.SUBCLASS_RESPONSIBILITY)
};
shindig.Gadget.prototype.getUserPrefsDialogContent=function(A){throw Error(shindig.errors.SUBCLASS_RESPONSIBILITY)
};
shindig.Gadget.prototype.getMainContent=function(A){throw Error(shindig.errors.SUBCLASS_RESPONSIBILITY)
};
shindig.Gadget.prototype.finishRender=function(A){throw Error(shindig.errors.SUBCLASS_RESPONSIBILITY)
};
shindig.Gadget.prototype.getAdditionalParams=function(){return""
};
shindig.BaseIfrGadget=function(A){shindig.Gadget.call(this,A);
this.serverBase_="/gadgets/";
this.queryIfrGadgetType_()
};
shindig.BaseIfrGadget.inherits(shindig.Gadget);
shindig.BaseIfrGadget.prototype.GADGET_IFRAME_PREFIX_="remote_iframe_";
shindig.BaseIfrGadget.prototype.CONTAINER="default";
shindig.BaseIfrGadget.prototype.cssClassGadget="gadgets-gadget";
shindig.BaseIfrGadget.prototype.cssClassTitleBar="gadgets-gadget-title-bar";
shindig.BaseIfrGadget.prototype.cssClassTitle="gadgets-gadget-title";
shindig.BaseIfrGadget.prototype.cssClassTitleButtonBar="gadgets-gadget-title-button-bar";
shindig.BaseIfrGadget.prototype.cssClassGadgetUserPrefsDialog="gadgets-gadget-user-prefs-dialog";
shindig.BaseIfrGadget.prototype.cssClassGadgetUserPrefsDialogActionBar="gadgets-gadget-user-prefs-dialog-action-bar";
shindig.BaseIfrGadget.prototype.cssClassTitleButton="gadgets-gadget-title-button";
shindig.BaseIfrGadget.prototype.cssClassGadgetContent="gadgets-gadget-content";
shindig.BaseIfrGadget.prototype.rpcToken=(2147483647*Math.random())|0;
shindig.BaseIfrGadget.prototype.rpcRelay="../container/rpc_relay.html";
shindig.BaseIfrGadget.prototype.getTitleBarContent=function(A){var B=this.hasViewablePrefs_()?'<a href="#" onclick="shindig.container.getGadget('+this.id+').handleOpenUserPrefsDialog();return false;" class="'+this.cssClassTitleButton+'">settings</a> ':"";
A('<div id="'+this.cssClassTitleBar+"-"+this.id+'" class="'+this.cssClassTitleBar+'"><span id="'+this.getIframeId()+'_title" class="'+this.cssClassTitle+'">'+(this.title?this.title:"Title")+'</span> | <span class="'+this.cssClassTitleButtonBar+'">'+B+'<a href="#" onclick="shindig.container.getGadget('+this.id+').handleToggle();return false;" class="'+this.cssClassTitleButton+'">toggle</a></span></div>')
};
shindig.BaseIfrGadget.prototype.getUserPrefsDialogContent=function(A){A('<div id="'+this.getUserPrefsDialogId()+'" class="'+this.cssClassGadgetUserPrefsDialog+'"></div>')
};
shindig.BaseIfrGadget.prototype.setServerBase=function(A){this.serverBase_=A
};
shindig.BaseIfrGadget.prototype.getServerBase=function(){return this.serverBase_
};
shindig.BaseIfrGadget.prototype.getMainContent=function(A){var B=this;
window.setTimeout(function(){B.getMainContent(A)
},0)
};
shindig.BaseIfrGadget.prototype.getIframeId=function(){return this.GADGET_IFRAME_PREFIX_+this.id
};
shindig.BaseIfrGadget.prototype.getUserPrefsDialogId=function(){return this.getIframeId()+"_userPrefsDialog"
};
shindig.BaseIfrGadget.prototype.getUserPrefsParams=function(){var B="";
for(var A in this.getUserPrefs()){B+="&up_"+encodeURIComponent(A)+"="+encodeURIComponent(this.getUserPrefValue(A))
}return B
};
shindig.BaseIfrGadget.prototype.handleToggle=function(){var B=document.getElementById(this.getIframeId());
if(B){var A=B.parentNode;
var C=A.style.display;
A.style.display=C?"":"none"
}};
shindig.BaseIfrGadget.prototype.hasViewablePrefs_=function(){for(var B in this.getUserPrefs()){var A=this.userPrefs[B];
if(A.type!="hidden"){return true
}}return false
};
shindig.BaseIfrGadget.prototype.handleOpenUserPrefsDialog=function(){if(this.userPrefsDialogContentLoaded){this.showUserPrefsDialog()
}else{var C=this;
var B="ig_callback_"+this.id;
window[B]=function(D){C.userPrefsDialogContentLoaded=true;
C.buildUserPrefsDialog(D);
C.showUserPrefsDialog()
};
var A=document.createElement("script");
A.src="http://www.gmodules.com/ig/gadgetsettings?mid="+this.id+"&output=js"+this.getUserPrefsParams()+"&url="+this.specUrl;
document.body.appendChild(A)
}};
shindig.BaseIfrGadget.prototype.buildUserPrefsDialog=function(A){var B=document.getElementById(this.getUserPrefsDialogId());
B.innerHTML=A+'<div class="'+this.cssClassGadgetUserPrefsDialogActionBar+'"><input type="button" value="Save" onclick="shindig.container.getGadget('+this.id+').handleSaveUserPrefs()"> <input type="button" value="Cancel" onclick="shindig.container.getGadget('+this.id+').handleCancelUserPrefs()"></div>';
B.childNodes[0].style.display=""
};
shindig.BaseIfrGadget.prototype.showUserPrefsDialog=function(A){var B=document.getElementById(this.getUserPrefsDialogId());
B.style.display=(A||A===undefined)?"":"none"
};
shindig.BaseIfrGadget.prototype.hideUserPrefsDialog=function(){this.showUserPrefsDialog(false)
};
shindig.BaseIfrGadget.prototype.handleSaveUserPrefs=function(){this.hideUserPrefsDialog();
var A=document.getElementById("m_"+this.id+"_numfields").value;
for(var D=0;
D<A;
D++){var B=document.getElementById("m_"+this.id+"_"+D);
var F="m_"+this.id+"_up_";
var C=B.name.substring(F.length);
var E=B.value;
this.userPrefs[C].value=E
}this.saveUserPrefs();
this.refresh()
};
shindig.BaseIfrGadget.prototype.handleCancelUserPrefs=function(){this.hideUserPrefsDialog()
};
shindig.BaseIfrGadget.prototype.refresh=function(){var A=this.getIframeId();
document.getElementById(A).src=this.getIframeUrl()
};
shindig.BaseIfrGadget.prototype.queryIfrGadgetType_=function(){var C={context:{country:"default",language:"default",view:"profile",container:"default"},gadgets:[{url:this.specUrl,moduleId:1}]};
var B={CONTENT_TYPE:"JSON",METHOD:"POST",POST_DATA:gadgets.json.stringify(C)};
var A=this.serverBase_+"metadata?st="+this.secureToken;
gadgets.io.makeNonProxiedRequest(A,D,B,"application/javascript");
var E=this;
function D(J){var K=false;
var F=J.data.gadgets[0].features;
if(F!=null){for(var H=0;
H<F.length;
H++){if(F[H]==="pubsub-2"){K=true;
break
}}}var I=K?shindig.OAAIfrGadget:shindig.IfrGadget;
for(var G in I){if(I.hasOwnProperty(G)){E[G]=I[G]
}}}};
shindig.IfrGadget={getMainContent:function(A){var B=this.getIframeId();
gadgets.rpc.setRelayUrl(B,this.serverBase_+this.rpcRelay);
gadgets.rpc.setAuthToken(B,this.rpcToken);
A('<div class="'+this.cssClassGadgetContent+'"><iframe id="'+B+'" name="'+B+'" class="'+this.cssClassGadget+'" src="about:blank" frameborder="no" scrolling="no"'+(this.height?' height="'+this.height+'"':"")+(this.width?' width="'+this.width+'"':"")+"></iframe></div>")
},finishRender:function(A){window.frames[this.getIframeId()].location=this.getIframeUrl()
},getIframeUrl:function(){return this.serverBase_+"ifr?container="+this.CONTAINER+"&mid="+this.id+"&nocache="+shindig.container.nocache_+"&country="+shindig.container.country_+"&lang="+shindig.container.language_+"&view="+shindig.container.view_+(this.specVersion?"&v="+this.specVersion:"")+(shindig.container.parentUrl_?"&parent="+encodeURIComponent(shindig.container.parentUrl_):"")+(this.debug?"&debug=1":"")+this.getAdditionalParams()+this.getUserPrefsParams()+(this.secureToken?"&st="+this.secureToken:"")+"&url="+encodeURIComponent(this.specUrl)+"#rpctoken="+this.rpcToken+(this.viewParams?"&view-params="+encodeURIComponent(gadgets.json.stringify(this.viewParams)):"")+(this.hashData?"&"+this.hashData:"")
}};
shindig.OAAIfrGadget={getMainContent:function(A){A('<div id="'+this.cssClassGadgetContent+"-"+this.id+'" class="'+this.cssClassGadgetContent+'"></div>')
},finishRender:function(A){var B={className:this.cssClassGadget,frameborder:"no",scrolling:"no"};
if(this.height){B.height=this.height
}if(this.width){B.width=this.width
}new OpenAjax.hub.IframeContainer(gadgets.pubsub2router.hub,this.getIframeId(),{Container:{onSecurityAlert:function(D,C){gadgets.error("Security error for container "+D.getClientID()+" : "+C);
D.getIframe().src="about:blank"
}},IframeContainer:{parent:document.getElementById(this.cssClassGadgetContent+"-"+this.id),uri:this.getIframeUrl(),tunnelURI:shindig.uri(this.serverBase_+this.rpcRelay).resolve(shindig.uri(window.location.href)),iframeAttrs:B}})
},getIframeUrl:function(){return this.serverBase_+"ifr?container="+this.CONTAINER+"&mid="+this.id+"&nocache="+shindig.container.nocache_+"&country="+shindig.container.country_+"&lang="+shindig.container.language_+"&view="+shindig.container.view_+(this.specVersion?"&v="+this.specVersion:"")+(this.debug?"&debug=1":"")+this.getAdditionalParams()+this.getUserPrefsParams()+(this.secureToken?"&st="+this.secureToken:"")+"&url="+encodeURIComponent(this.specUrl)+(this.viewParams?"&view-params="+encodeURIComponent(gadgets.json.stringify(this.viewParams)):"")+(this.hashData?"#"+this.hashData:"")
}};
shindig.Container=function(){this.gadgets_={};
this.parentUrl_=document.location.href+"://"+document.location.host;
this.country_="ALL";
this.language_="ALL";
this.view_="profile";
this.nocache_=1;
this.maxheight_=2147483647
};
shindig.Container.inherits(shindig.Extensible);
shindig.Container.prototype.gadgetClass=shindig.Gadget;
shindig.Container.prototype.userPrefStore=new shindig.DefaultUserPrefStore();
shindig.Container.prototype.gadgetService=new shindig.GadgetService();
shindig.Container.prototype.layoutManager=new shindig.StaticLayoutManager();
shindig.Container.prototype.setParentUrl=function(A){this.parentUrl_=A
};
shindig.Container.prototype.setCountry=function(A){this.country_=A
};
shindig.Container.prototype.setNoCache=function(A){this.nocache_=A
};
shindig.Container.prototype.setLanguage=function(A){this.language_=A
};
shindig.Container.prototype.setView=function(A){this.view_=A
};
shindig.Container.prototype.setMaxHeight=function(A){this.maxheight_=A
};
shindig.Container.prototype.getGadgetKey_=function(A){return"gadget_"+A
};
shindig.Container.prototype.getGadget=function(A){return this.gadgets_[this.getGadgetKey_(A)]
};
shindig.Container.prototype.createGadget=function(A){return new this.gadgetClass(A)
};
shindig.Container.prototype.addGadget=function(A){A.id=this.getNextGadgetInstanceId();
this.gadgets_[this.getGadgetKey_(A.id)]=A
};
shindig.Container.prototype.addGadgets=function(A){for(var B=0;
B<A.length;
B++){this.addGadget(A[B])
}};
shindig.Container.prototype.renderGadgets=function(){for(var A in this.gadgets_){this.renderGadget(this.gadgets_[A])
}};
shindig.Container.prototype.renderGadget=function(A){throw Error(shindig.errors.SUBCLASS_RESPONSIBILITY)
};
shindig.Container.prototype.nextGadgetInstanceId_=0;
shindig.Container.prototype.getNextGadgetInstanceId=function(){return this.nextGadgetInstanceId_++
};
shindig.Container.prototype.refreshGadgets=function(){for(var A in this.gadgets_){this.gadgets_[A].refresh()
}};
shindig.IfrContainer=function(){shindig.Container.call(this)
};
shindig.IfrContainer.inherits(shindig.Container);
shindig.IfrContainer.prototype.gadgetClass=shindig.BaseIfrGadget;
shindig.IfrContainer.prototype.gadgetService=new shindig.IfrGadgetService();
shindig.IfrContainer.prototype.setParentUrl=function(A){if(!A.match(/^http[s]?:\/\//)){A=document.location.href.match(/^[^?#]+\//)[0]+A
}this.parentUrl_=A
};
shindig.IfrContainer.prototype.renderGadget=function(B){var A=this.layoutManager.getGadgetChrome(B);
B.render(A)
};
shindig.container=new shindig.IfrContainer();;
if(gadgets&&gadgets.rpc){osapi._handleGadgetRpcMethod=function(A){var F=new Array(A.length);
var E=0;
var H=this.callback;
var B=function(K,J){J({})
};
for(var D=0;
D<A.length;
D++){var G=osapi;
if(A[D].method.indexOf("_")==-1){var I=A[D].method.split(".");
for(var C=0;
C<I.length;
C++){if(G.hasOwnProperty(I[C])){G=G[I[C]]
}else{G=B;
break
}}}else{G=B
}G(A[D].params,function(J){return function(K){F[J]={id:A[J].id,data:K};
E++;
if(E==A.length){H(F)
}}
}(D))
}};
osapi.container={};
osapi.container.listMethods=function(A,C){var B=[];
recurseNames(osapi,"",5,B);
C(B)
};
function recurseNames(C,D,E,B){if(E==0){return 
}for(var F in C){if(C.hasOwnProperty(F)){if(F.indexOf("_")==-1){var A=typeof (C[F]);
if(A=="function"){B.push(D+F)
}else{if(A=="object"){recurseNames(C[F],D+F+".",E-1,B)
}}}}}}gadgets.rpc.register("osapi._handleGadgetRpcMethod",osapi._handleGadgetRpcMethod)
};;
var OpenAjax=OpenAjax||{};
if(!OpenAjax.hub){OpenAjax.hub=function(){var B={};
var A="org.openajax.hub.";
return{implementer:"http://openajax.org",implVersion:"2.0.4",specVersion:"2.0",implExtraData:{},libraries:B,registerLibrary:function(F,E,D,C){B[F]={prefix:F,namespaceURI:E,version:D,extraData:C};
this.publish(A+"registerLibrary",B[F])
},unregisterLibrary:function(C){this.publish(A+"unregisterLibrary",B[C]);
delete B[C]
}}
}();
OpenAjax.hub.Error={BadParameters:"OpenAjax.hub.Error.BadParameters",Disconnected:"OpenAjax.hub.Error.Disconnected",Duplicate:"OpenAjax.hub.Error.Duplicate",NoContainer:"OpenAjax.hub.Error.NoContainer",NoSubscription:"OpenAjax.hub.Error.NoSubscription",NotAllowed:"OpenAjax.hub.Error.NotAllowed",WrongProtocol:"OpenAjax.hub.Error.WrongProtocol",IncompatBrowser:"OpenAjax.hub.Error.IncompatBrowser"};
OpenAjax.hub.SecurityAlert={LoadTimeout:"OpenAjax.hub.SecurityAlert.LoadTimeout",FramePhish:"OpenAjax.hub.SecurityAlert.FramePhish",ForgedMsg:"OpenAjax.hub.SecurityAlert.ForgedMsg"};
OpenAjax.hub._debugger=function(){};
OpenAjax.hub.ManagedHub=function(B){if(!B||!B.onPublish||!B.onSubscribe){throw new Error(OpenAjax.hub.Error.BadParameters)
}this._p=B;
this._onUnsubscribe=B.onUnsubscribe?B.onUnsubscribe:null;
this._scope=B.scope||window;
if(B.log){var A=this;
this._log=function(D){try{B.log.call(A._scope,"ManagedHub: "+D)
}catch(C){OpenAjax.hub._debugger()
}}
}else{this._log=function(){}
}this._subscriptions={c:{},s:null};
this._containers={};
this._seq=0;
this._active=true;
this._isPublishing=false;
this._pubQ=[]
};
OpenAjax.hub.ManagedHub.prototype.subscribeForClient=function(A,B,C){this._assertConn();
if(this._invokeOnSubscribe(B,A)){return this._subscribe(B,this._sendToClient,this,{c:A,sid:C})
}throw new Error(OpenAjax.hub.Error.NotAllowed)
};
OpenAjax.hub.ManagedHub.prototype.unsubscribeForClient=function(A,B){this._unsubscribe(B);
this._invokeOnUnsubscribe(A,B)
};
OpenAjax.hub.ManagedHub.prototype.publishForClient=function(A,B,C){this._assertConn();
this._publish(B,C,A)
};
OpenAjax.hub.ManagedHub.prototype.disconnect=function(){this._active=false;
for(var A in this._containers){this.removeContainer(this._containers[A])
}};
OpenAjax.hub.ManagedHub.prototype.getContainer=function(B){var A=this._containers[B];
return A?A:null
};
OpenAjax.hub.ManagedHub.prototype.listContainers=function(){var A=[];
for(var B in this._containers){A.push(this._containers[B])
}return A
};
OpenAjax.hub.ManagedHub.prototype.addContainer=function(B){this._assertConn();
var A=B.getClientID();
if(this._containers[A]){throw new Error(OpenAjax.hub.Error.Duplicate)
}this._containers[A]=B
};
OpenAjax.hub.ManagedHub.prototype.removeContainer=function(B){var A=B.getClientID();
if(!this._containers[A]){throw new Error(OpenAjax.hub.Error.NoContainer)
}B.remove();
delete this._containers[A]
};
OpenAjax.hub.ManagedHub.prototype.subscribe=function(B,E,D,H,C){this._assertConn();
this._assertSubTopic(B);
if(!E){throw new Error(OpenAjax.hub.Error.BadParameters)
}D=D||window;
if(!this._invokeOnSubscribe(B,null)){this._invokeOnComplete(H,D,null,false,OpenAjax.hub.Error.NotAllowed);
return 
}var F=this;
function G(J,K,M,I){if(F._invokeOnPublish(J,K,I,null)){try{E.call(D,J,K,C)
}catch(L){OpenAjax.hub._debugger();
F._log("caught error from onData callback to Hub.subscribe(): "+L.message)
}}}var A=this._subscribe(B,G,D,C);
this._invokeOnComplete(H,D,A,true);
return A
};
OpenAjax.hub.ManagedHub.prototype.publish=function(A,B){this._assertConn();
this._assertPubTopic(A);
this._publish(A,B,null)
};
OpenAjax.hub.ManagedHub.prototype.unsubscribe=function(A,C,B){this._assertConn();
if(!A){throw new Error(OpenAjax.hub.Error.BadParameters)
}this._unsubscribe(A);
this._invokeOnUnsubscribe(null,A);
this._invokeOnComplete(C,B,A,true)
};
OpenAjax.hub.ManagedHub.prototype.isConnected=function(){return this._active
};
OpenAjax.hub.ManagedHub.prototype.getScope=function(){return this._scope
};
OpenAjax.hub.ManagedHub.prototype.getSubscriberData=function(C){this._assertConn();
var D=C.split(".");
var A=D.pop();
var B=this._getSubscriptionObject(this._subscriptions,D,0,A);
if(B){return B.data
}throw new Error(OpenAjax.hub.Error.NoSubscription)
};
OpenAjax.hub.ManagedHub.prototype.getSubscriberScope=function(C){this._assertConn();
var D=C.split(".");
var A=D.pop();
var B=this._getSubscriptionObject(this._subscriptions,D,0,A);
if(B){return B.scope
}throw new Error(OpenAjax.hub.Error.NoSubscription)
};
OpenAjax.hub.ManagedHub.prototype.getParameters=function(){return this._p
};
OpenAjax.hub.ManagedHub.prototype._sendToClient=function(B,C,D,A){if(!this.isConnected()){return 
}if(this._invokeOnPublish(B,C,A,D.c)){D.c.sendToClient(B,C,D.sid)
}};
OpenAjax.hub.ManagedHub.prototype._assertConn=function(){if(!this.isConnected()){throw new Error(OpenAjax.hub.Error.Disconnected)
}};
OpenAjax.hub.ManagedHub.prototype._assertPubTopic=function(A){if(!A||A===""||(A.indexOf("*")!=-1)||(A.indexOf("..")!=-1)||(A.charAt(0)==".")||(A.charAt(A.length-1)==".")){throw new Error(OpenAjax.hub.Error.BadParameters)
}};
OpenAjax.hub.ManagedHub.prototype._assertSubTopic=function(B){if(!B){throw new Error(OpenAjax.hub.Error.BadParameters)
}var E=B.split(".");
var A=E.length;
for(var C=0;
C<A;
C++){var D=E[C];
if((D==="")||((D.indexOf("*")!=-1)&&(D!="*")&&(D!="**"))){throw new Error(OpenAjax.hub.Error.BadParameters)
}if((D=="**")&&(C<A-1)){throw new Error(OpenAjax.hub.Error.BadParameters)
}}};
OpenAjax.hub.ManagedHub.prototype._invokeOnComplete=function(C,A,B,F,E){if(C){try{A=A||window;
C.call(A,B,F,E)
}catch(D){OpenAjax.hub._debugger();
this._log("caught error from onComplete callback: "+D.message)
}}};
OpenAjax.hub.ManagedHub.prototype._invokeOnPublish=function(B,C,A,E){try{return this._p.onPublish.call(this._scope,B,C,A,E)
}catch(D){OpenAjax.hub._debugger();
this._log("caught error from onPublish callback to constructor: "+D.message)
}return false
};
OpenAjax.hub.ManagedHub.prototype._invokeOnSubscribe=function(B,A){try{return this._p.onSubscribe.call(this._scope,B,A)
}catch(C){OpenAjax.hub._debugger();
this._log("caught error from onSubscribe callback to constructor: "+C.message)
}return false
};
OpenAjax.hub.ManagedHub.prototype._invokeOnUnsubscribe=function(A,C){if(this._onUnsubscribe){var B=C.slice(0,C.lastIndexOf("."));
try{this._onUnsubscribe.call(this._scope,B,A)
}catch(D){OpenAjax.hub._debugger();
this._log("caught error from onUnsubscribe callback to constructor: "+D.message)
}}};
OpenAjax.hub.ManagedHub.prototype._subscribe=function(A,E,D,C){var F=A+"."+this._seq;
var B={scope:D,cb:E,data:C,sid:this._seq++};
var G=A.split(".");
this._recursiveSubscribe(this._subscriptions,G,0,B);
return F
};
OpenAjax.hub.ManagedHub.prototype._recursiveSubscribe=function(A,E,B,D){var C=E[B];
if(B==E.length){D.next=A.s;
A.s=D
}else{if(typeof A.c=="undefined"){A.c={}
}if(typeof A.c[C]=="undefined"){A.c[C]={c:{},s:null};
this._recursiveSubscribe(A.c[C],E,B+1,D)
}else{this._recursiveSubscribe(A.c[C],E,B+1,D)
}}};
OpenAjax.hub.ManagedHub.prototype._publish=function(B,D,A){if(this._isPublishing){this._pubQ.push({t:B,d:D,p:A});
return 
}this._safePublish(B,D,A);
while(this._pubQ.length>0){var C=this._pubQ.shift();
this._safePublish(C.t,C.d,C.p)
}};
OpenAjax.hub.ManagedHub.prototype._safePublish=function(B,C,A){this._isPublishing=true;
var D=B.split(".");
this._recursivePublish(this._subscriptions,D,0,B,C,A);
this._isPublishing=false
};
OpenAjax.hub.ManagedHub.prototype._recursivePublish=function(K,J,G,B,C,E){if(typeof K!="undefined"){var D;
if(G==J.length){D=K
}else{this._recursivePublish(K.c[J[G]],J,G+1,B,C,E);
this._recursivePublish(K.c["*"],J,G+1,B,C,E);
D=K.c["**"]
}if(typeof D!="undefined"){var A=D.s;
while(A){var I=A.scope;
var F=A.cb;
var H=A.data;
if(typeof F=="string"){F=I[F]
}F.call(I,B,C,H,E);
A=A.next
}}}};
OpenAjax.hub.ManagedHub.prototype._unsubscribe=function(B){var C=B.split(".");
var A=C.pop();
if(!this._recursiveUnsubscribe(this._subscriptions,C,0,A)){throw new Error(OpenAjax.hub.Error.NoSubscription)
}};
OpenAjax.hub.ManagedHub.prototype._recursiveUnsubscribe=function(I,H,E,D){if(typeof I=="undefined"){return false
}if(E<H.length){var C=I.c[H[E]];
if(!C){return false
}this._recursiveUnsubscribe(C,H,E+1,D);
if(!C.s){for(var F in C.c){return true
}delete I.c[H[E]]
}}else{var B=I.s;
var A=null;
var G=false;
while(B){if(D==B.sid){G=true;
if(B==I.s){I.s=B.next
}else{A.next=B.next
}break
}A=B;
B=B.next
}if(!G){return false
}}return true
};
OpenAjax.hub.ManagedHub.prototype._getSubscriptionObject=function(A,F,C,B){if(typeof A!="undefined"){if(C<F.length){var D=A.c[F[C]];
return this._getSubscriptionObject(D,F,C+1,B)
}var E=A.s;
while(E){if(B==E.sid){return E
}E=E.next
}}return null
};
OpenAjax.hub._hub=new OpenAjax.hub.ManagedHub({onSubscribe:function(A,B){return true
},onPublish:function(B,C,A,D){return true
}});
OpenAjax.hub.subscribe=function(A,D,C,B){if(typeof D==="string"){C=C||window;
D=C[D]||null
}return OpenAjax.hub._hub.subscribe(A,D,C,null,B)
};
OpenAjax.hub.unsubscribe=function(A){return OpenAjax.hub._hub.unsubscribe(A)
};
OpenAjax.hub.publish=function(A,B){OpenAjax.hub._hub.publish(A,B)
};
OpenAjax.hub.registerLibrary("OpenAjax","http://openajax.org/hub","2.0",{})
};;
var OpenAjax=OpenAjax||{};
OpenAjax.hub=OpenAjax.hub||{};
OpenAjax.gadgets=typeof OpenAjax.gadgets==="object"?OpenAjax.gadgets:typeof gadgets==="object"?gadgets:{};
OpenAjax.gadgets.rpctx=OpenAjax.gadgets.rpctx||{};
(function(){if(typeof gadgets==="undefined"){if(typeof oaaConfig==="undefined"){var scripts=document.getElementsByTagName("script");
var reHub=/openajax(?:managedhub-(?:all|core).*|-mashup)\.js$/i;
for(var i=scripts.length-1;
i>=0;
i--){var src=scripts[i].getAttribute("src");
if(!src){continue
}var m=src.match(reHub);
if(m){var config=scripts[i].getAttribute("oaaConfig");
if(config){try{oaaConfig=eval("({ "+config+" })")
}catch(e){}}break
}}}if(typeof oaaConfig!=="undefined"&&oaaConfig.gadgetsGlobal){gadgets=OpenAjax.gadgets
}}})();
if(!OpenAjax.hub.IframeContainer){(function(){OpenAjax.hub.IframeContainer=function(I,O,S){E(arguments);
var M=this;
var D=S.Container.scope||window;
var T=false;
var R={};
var N;
var P;
var K=S.IframeContainer.timeout||15000;
var C;
if(S.Container.log){var J=function(V){try{S.Container.log.call(D,"IframeContainer::"+O+": "+V)
}catch(U){OpenAjax.hub._debugger()
}}
}else{J=function(){}
}this._init=function(){I.addContainer(this);
P=OpenAjax.hub.IframeContainer._rpcRouter.add(O,this);
N=A(S,D,J);
var V=null;
var U=OpenAjax.gadgets.rpc.getRelayChannel();
if(S.IframeContainer.tunnelURI){if(U!=="wpm"&&U!=="nix"){throw new Error(OpenAjax.hub.Error.IncompatBrowser)
}}else{J("WARNING: Parameter 'IframeContaienr.tunnelURI' not specified. Connection will not be fully secure.");
if(U==="rmr"){V=OpenAjax.gadgets.rpc.getOrigin(S.IframeContainer.uri)+"/robots.txt"
}}F();
OpenAjax.gadgets.rpc.setupReceiver(P,V);
L()
};
this.sendToClient=function(V,W,U){OpenAjax.gadgets.rpc.call(P,"openajax.pubsub",null,"pub",V,W,U)
};
this.remove=function(){H();
clearTimeout(C);
OpenAjax.gadgets.rpc.removeReceiver(P);
var U=document.getElementById(P);
U.parentNode.removeChild(U);
OpenAjax.hub.IframeContainer._rpcRouter.remove(P)
};
this.isConnected=function(){return T
};
this.getClientID=function(){return O
};
this.getPartnerOrigin=function(){if(T){var U=OpenAjax.gadgets.rpc.getReceiverOrigin(P);
if(U){return(/^([a-zA-Z]+:\/\/[^:]+).*/.exec(U)[1])
}}return null
};
this.getParameters=function(){return S
};
this.getHub=function(){return I
};
this.getIframe=function(){return document.getElementById(P)
};
function E(U){var V=U[0],X=U[1],W=U[2];
if(!V||!X||!W||!W.Container||!W.Container.onSecurityAlert||!W.IframeContainer||!W.IframeContainer.parent||!W.IframeContainer.uri){throw new Error(OpenAjax.hub.Error.BadParameters)
}}this._handleIncomingRPC=function(Z,V,X){switch(Z){case"pub":I.publishForClient(M,V,X);
break;
case"sub":var U="";
try{R[X]=I.subscribeForClient(M,V,X)
}catch(Y){U=Y.message
}return U;
case"uns":var W=R[X];
I.unsubscribeForClient(M,W);
delete R[X];
return X;
case"con":G();
return true;
case"dis":L();
H();
if(S.Container.onDisconnect){try{S.Container.onDisconnect.call(D,M)
}catch(Y){OpenAjax.hub._debugger();
J("caught error from onDisconnect callback to constructor: "+Y.message)
}}return true
}};
this._onSecurityAlert=function(U){Q(B[U])
};
function F(){var Y=document.createElement("span");
S.IframeContainer.parent.appendChild(Y);
var W='<iframe id="'+P+'" name="'+P+'" src="javascript:\'<html></html>\'"';
var a="";
var V=S.IframeContainer.iframeAttrs;
if(V){for(var U in V){switch(U){case"style":for(var X in V.style){a+=X+":"+V.style[X]+";"
}break;
case"className":W+=' class="'+V[U]+'"';
break;
default:W+=" "+U+'="'+V[U]+'"'
}}}a+="visibility:hidden;";
W+=' style="'+a+'"></iframe>';
Y.innerHTML=W;
var Z=S.IframeContainer.tunnelURI;
document.getElementById(P).src=S.IframeContainer.uri+"#rpctoken="+N+(Z?"&parent="+encodeURIComponent(Z)+"&forcesecure=true":"&oaaParent="+encodeURIComponent(OpenAjax.gadgets.rpc.getOrigin(window.location.href)))
}function G(){function U(V){if(V){T=true;
clearTimeout(C);
document.getElementById(P).style.visibility="visible";
if(S.Container.onConnect){try{S.Container.onConnect.call(D,M)
}catch(W){OpenAjax.hub._debugger();
J("caught error from onConnect callback to constructor: "+W.message)
}}}}OpenAjax.gadgets.rpc.call(P,"openajax.pubsub",U,"cmd","con")
}function H(){if(T){T=false;
document.getElementById(P).style.visibility="hidden";
for(var U in R){I.unsubscribeForClient(M,R[U])
}R={}
}}function Q(U){try{S.Container.onSecurityAlert.call(D,M,U)
}catch(V){OpenAjax.hub._debugger();
J("caught error from onSecurityAlert callback to constructor: "+V.message)
}}function L(){C=setTimeout(function(){Q(OpenAjax.hub.SecurityAlert.LoadTimeout);
M._handleIncomingRPC=function(){}
},K)
}this._init()
};
OpenAjax.hub.IframeHubClient=function(F){if(!F||!F.HubClient||!F.HubClient.onSecurityAlert){throw new Error(OpenAjax.hub.Error.BadParameters)
}var E=this;
var M=F.HubClient.scope||window;
var I=false;
var G={};
var H=0;
var K;
if(F.HubClient.log){var J=function(O){try{F.HubClient.log.call(M,"IframeHubClient::"+K+": "+O)
}catch(N){OpenAjax.hub._debugger()
}}
}else{J=function(){}
}this._init=function(){var O=OpenAjax.gadgets.util.getUrlParameters();
if(!O.parent){var P=O.oaaParent+"/robots.txt";
OpenAjax.gadgets.rpc.setupReceiver("..",P)
}if(F.IframeHubClient&&F.IframeHubClient.requireParentVerifiable&&OpenAjax.gadgets.rpc.getReceiverOrigin("..")===null){OpenAjax.gadgets.rpc.removeReceiver("..");
throw new Error(OpenAjax.hub.Error.IncompatBrowser)
}OpenAjax.hub.IframeContainer._rpcRouter.add("..",this);
var N=OpenAjax.gadgets.rpc.RPC_ID;
if(!N){throw new Error(OpenAjax.hub.Error.WrongProtocol)
}K=N.substr(N.indexOf("_")+1)
};
this.connect=function(O,N){if(I){throw new Error(OpenAjax.hub.Error.Duplicate)
}function P(Q){if(Q){I=true;
if(O){try{O.call(N||window,E,true)
}catch(R){OpenAjax.hub._debugger();
J("caught error from onComplete callback to connect(): "+R.message)
}}}}OpenAjax.gadgets.rpc.call("..","openajax.pubsub",P,"con")
};
this.disconnect=function(O,N){if(!I){throw new Error(OpenAjax.hub.Error.Disconnected)
}I=false;
var P=null;
if(O){P=function(Q){try{O.call(N||window,E,true)
}catch(R){OpenAjax.hub._debugger();
J("caught error from onComplete callback to disconnect(): "+R.message)
}}
}OpenAjax.gadgets.rpc.call("..","openajax.pubsub",P,"dis")
};
this.getPartnerOrigin=function(){if(I){var N=OpenAjax.gadgets.rpc.getReceiverOrigin("..");
if(N){return(/^([a-zA-Z]+:\/\/[^:]+).*/.exec(N)[1])
}}return null
};
this.getClientID=function(){return K
};
this.subscribe=function(O,R,Q,S,P){D();
C(O);
if(!R){throw new Error(OpenAjax.hub.Error.BadParameters)
}Q=Q||window;
var N=""+H++;
G[N]={cb:R,sc:Q,d:P};
function T(U){if(U!==""){delete G[N]
}if(S){try{S.call(Q,N,U==="",U)
}catch(V){OpenAjax.hub._debugger();
J("caught error from onComplete callback to subscribe(): "+V.message)
}}}OpenAjax.gadgets.rpc.call("..","openajax.pubsub",T,"sub",O,N);
return N
};
this.publish=function(N,O){D();
L(N);
OpenAjax.gadgets.rpc.call("..","openajax.pubsub",null,"pub",N,O)
};
this.unsubscribe=function(N,P,O){D();
if(!N){throw new Error(OpenAjax.hub.Error.BadParameters)
}if(!G[N]||G[N].uns){throw new Error(OpenAjax.hub.Error.NoSubscription)
}G[N].uns=true;
function Q(R){delete G[N];
if(P){try{P.call(O||window,N,true)
}catch(S){OpenAjax.hub._debugger();
J("caught error from onComplete callback to unsubscribe(): "+S.message)
}}}OpenAjax.gadgets.rpc.call("..","openajax.pubsub",Q,"uns",null,N)
};
this.isConnected=function(){return I
};
this.getScope=function(){return M
};
this.getSubscriberData=function(N){D();
if(G[N]){return G[N].d
}throw new Error(OpenAjax.hub.Error.NoSubscription)
};
this.getSubscriberScope=function(N){D();
if(G[N]){return G[N].sc
}throw new Error(OpenAjax.hub.Error.NoSubscription)
};
this.getParameters=function(){return F
};
this._handleIncomingRPC=function(R,O,P,N){if(R==="pub"){if(G[N]&&!G[N].uns){try{G[N].cb.call(G[N].sc,O,P,G[N].d)
}catch(Q){OpenAjax.hub._debugger();
J("caught error from onData callback to subscribe(): "+Q.message)
}}}if(O==="con"){return true
}return false
};
function D(){if(!I){throw new Error(OpenAjax.hub.Error.Disconnected)
}}function C(O){if(!O){throw new Error(OpenAjax.hub.Error.BadParameters)
}var R=O.split(".");
var N=R.length;
for(var P=0;
P<N;
P++){var Q=R[P];
if((Q==="")||((Q.indexOf("*")!=-1)&&(Q!="*")&&(Q!="**"))){throw new Error(OpenAjax.hub.Error.BadParameters)
}if((Q=="**")&&(P<N-1)){throw new Error(OpenAjax.hub.Error.BadParameters)
}}}function L(N){if(!N||N===""||(N.indexOf("*")!=-1)||(N.indexOf("..")!=-1)||(N.charAt(0)==".")||(N.charAt(N.length-1)==".")){throw new Error(OpenAjax.hub.Error.BadParameters)
}}this._init()
};
OpenAjax.hub.IframeContainer._rpcRouter=function(){var E={};
function D(){var F=E[this.f];
if(F){return F._handleIncomingRPC.apply(F,arguments)
}}function C(H,F){var G=E[H];
if(G){G._onSecurityAlert.call(G,F)
}}return{add:function(H,F){function G(J,I){if(J===".."){if(!E[".."]){E[".."]=I
}return 
}do{newID=((32767*Math.random())|0).toString(16)+"_"+J
}while(E[newID]);
E[newID]=I;
return newID
}OpenAjax.gadgets.rpc.register("openajax.pubsub",D);
OpenAjax.gadgets.rpc.config({securityCallback:C});
B[OpenAjax.gadgets.rpc.SEC_ERROR_LOAD_TIMEOUT]=OpenAjax.hub.SecurityAlert.LoadTimeout;
B[OpenAjax.gadgets.rpc.SEC_ERROR_FRAME_PHISH]=OpenAjax.hub.SecurityAlert.FramePhish;
B[OpenAjax.gadgets.rpc.SEC_ERROR_FORGED_MSG]=OpenAjax.hub.SecurityAlert.ForgedMsg;
this.add=G;
return G(H,F)
},remove:function(F){delete E[F]
}}
}();
var B={};
function A(H,E,D){if(!OpenAjax.hub.IframeContainer._prng){var C=new Date().getTime()+Math.random()+document.cookie;
OpenAjax.hub.IframeContainer._prng=OpenAjax._smash.crypto.newPRNG(C)
}var G=H.IframeContainer||H.IframeHubClient;
if(G&&G.seed){try{var J=G.seed.call(E);
OpenAjax.hub.IframeContainer._prng.addSeed(J)
}catch(F){OpenAjax.hub._debugger();
D("caught error from 'seed' callback: "+F.message)
}}var I=(G&&G.tokenLength)||6;
return OpenAjax.hub.IframeContainer._prng.nextRandomB64Str(I)
}})()
};;
if(typeof OpenAjax._smash=="undefined"){OpenAjax._smash={}
}OpenAjax._smash.crypto={strToWA:function(D,E){var C=Array();
var A=(1<<E)-1;
for(var B=0;
B<D.length*E;
B+=E){C[B>>5]|=(D.charCodeAt(B/E)&A)<<(32-E-B%32)
}return C
},hmac_sha1:function(D,G,F){var A=Array(16),C=Array(16);
for(var B=0;
B<16;
B++){A[B]=D[B]^909522486;
C[B]=D[B]^1549556828
}var E=this.sha1(A.concat(this.strToWA(G,F)),512+G.length*F);
return this.sha1(C.concat(E),512+160)
},newPRNG:function(A){var G=this;
if((typeof A!="string")||(A.length<12)){alert("WARNING: Seed length too short ...")
}var C=[43417,15926,18182,33130,9585,30800,49772,40144,47678,55453,4659,38181,65340,6787,54417,65301];
var B=[];
var E=0;
function F(H){return G.hmac_sha1(C,H,8)
}function D(H){var J=F(H);
for(var I=0;
I<5;
I++){B[I]^=J[I]
}}D(A);
return{addSeed:function(H){D(H)
},nextRandomOctets:function(H){var I=[];
while(H>0){E+=1;
var J=G.hmac_sha1(B,(E).toString(16),8);
for(i=0;
(i<20)&(H>0);
i++,H--){I.push((J[i>>2]>>(i%4))%256)
}}return I
},nextRandomB64Str:function(H){var L="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_";
var K=this.nextRandomOctets(H);
var I="";
for(var J=0;
J<H;
J++){I+=L.charAt(K[J]&63)
}return I
}}
},sha1:function(){var D=function(E,H){var G=(E&65535)+(H&65535);
var F=(E>>16)+(H>>16)+(G>>16);
return(F<<16)|(G&65535)
};
var C=function(E,F){return(E<<F)|(E>>>(32-F))
};
function B(F,E,H,G){if(F<20){return(E&H)|((~E)&G)
}if(F<40){return E^H^G
}if(F<60){return(E&H)|(E&G)|(H&G)
}return E^H^G
}function A(E){return(E<20)?1518500249:(E<40)?1859775393:(E<60)?-1894007588:-899497514
}return function(U,F){U[F>>5]|=128<<(24-F%32);
U[((F+64>>9)<<4)+15]=F;
var E=Array(80);
var S=1732584193;
var R=-271733879;
var Q=-1732584194;
var O=271733878;
var M=-1009589776;
for(var I=0;
I<U.length;
I+=16){var P=S;
var N=R;
var L=Q;
var K=O;
var J=M;
for(var H=0;
H<80;
H++){E[H]=((H<16)?U[I+H]:C(E[H-3]^E[H-8]^E[H-14]^E[H-16],1));
var G=D(D(C(P,5),B(H,N,L,K)),D(D(J,E[H]),A(H)));
J=K;
K=L;
L=C(N,30);
N=P;
P=G
}S=D(P,S);
R=D(N,R);
Q=D(L,Q);
O=D(K,O);
M=D(J,M)
}return Array(S,R,Q,O,M)
}
}()};;
gadgets.pubsub2router=function(){return{init:function(A){if(A.hub){this.hub=A.hub
}else{this.hub=new OpenAjax.hub.ManagedHub({onPublish:A.onPublish,onSubscribe:A.onSubscribe,onUnsubscribe:A.onUnsubscribe})
}}}
}();;
gadgets.config.init({"shindig.auth":{"authToken":"-1:-1:*::*:0:default"},"osapi":{"endPoints":["http://%host%/rpc"]},"osapi.services":{"gadgets.rpc":["container.listMethods"],"http://%host%/rpc":["samplecontainer.update","albums.supportedFields","albums.update","activities.delete","gadgets.metadata","activities.update","activities.supportedFields","mediaItems.create","albums.get","activities.get","http.put","activitystreams.create","messages.modify","appdata.get","messages.get","system.listMethods","samplecontainer.get","cache.invalidate","people.supportedFields","http.head","http.delete","messages.create","people.get","activitystreams.get","mediaItems.supportedFields","mediaItems.delete","albums.delete","activitystreams.update","mediaItems.update","messages.delete","appdata.update","gadgets.tokenSupportedFields","http.post","activities.create","samplecontainer.create","http.get","albums.create","appdata.delete","gadgets.token","appdata.create","activitystreams.delete","gadgets.supportedFields","mediaItems.get","activitystreams.supportedFields"]},"rpc":{"passReferrer":"c2p:query","parentRelayUrl":"/container/rpc_relay.html","useLegacyProtocol":false,"commSwf":"http://idc311-sciverse-shindig.elsevier.com/xpc.swf"},"core.io":{"proxyUrl":"//%host%/gadgets/proxy?container=default&refresh=%refresh%&url=%url%%rewriteMime%","jsonProxyUrl":"//%host%/gadgets/makeRequest"}});
