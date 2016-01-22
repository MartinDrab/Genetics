/***   Cytoscape code START  ***/
var cyjsFileName;
var jsonFileName;
var validCytoscapeNetwork;
var cyjsFileContent;
var jsonFileContent;
var cyElement;
var zoomSlider;
var isNodeEdge = true;
var currentElement;
var isFullScreen = false;
var environment = 0; //Environment Key values: 0-PROD, 1-CERT,2-DEV, 3-LOCAL
var queryStrings;
var currentIndex = -1;
var appName = "CytoscapeWeb";
var doi = ''; //doi=10.1016/S9999-9994(15)20790-8
var serverURL = [
'http://elsevier-apps.sciverse.com/',
'http://cert-elsevier-apps.sciverse.com/',
'http://dvc77761:26652/',
'http://138.12.195.20:8080/'
];

var isAppLoaded = false;

$(function(){ // on dom ready
	
	//STEP #5 - initialize Site catalyst impression on gadget makeMeVisible
	scPlatform = "SD";
	scGadgetName = "CI Cytoscape";
	
	queryStrings = getQueryVariable();	
	if(queryStrings){
		doi = (typeof queryStrings.doi === "undefined" || queryStrings.doi == '')?'':queryStrings.doi;
		environment = (typeof queryStrings.domain === "undefined" || queryStrings.domain == "")?environment:queryStrings.domain;
		currentIndex = (typeof queryStrings.index === "undefined" || queryStrings.index == "")?currentIndex:Number(queryStrings.index);
		isFullScreen = (typeof queryStrings.isFullScreen === "undefined" || queryStrings.isFullScreen == '')?isFullScreen:queryStrings.isFullScreen;	
	}	
	
	if(doi != ''){
		if(isFullScreen){
			$(".cyContainer").removeClass('cyGadgetView');
			$(".cyContainer").addClass('cyFullView');
			$("#appContainerId").addClass('appContainer');
			$("#headerDiv").css("display","none");
			$(".expIcon").css("display","none");
			$(".networkDiv").css('padding','8px 0 18px 0');
			$(".titlebar").css("display","block");
			$(".networkDiv").width("45%");
			$('html, body').css('overflow', 'visible');
			$(".nodeTxt").text("Show node labels");
			$(".edgeTxt").text("Show edge labels");
			$(".footerText").css("padding-left","75%");
			$("#invalidNetwork").css("width","auto");
			/*$("#invalidNetwork").css("padding-top","170px");*/
		}
		
		var fileUrl = serverURL[environment] + appName +"/CytoscapeFileServlet?doi="+ doi+"&timestamp="+getTimeStamp();
		
		//Start ajax request to fetch the cytscape network data from the BCR location
		$.ajax({
			url: fileUrl,		
			//force to handle it as text
			dataType: "json",
			success: function(data) {
				if(data && data.whiteListPII == true){
					console.log("data.whiteListPII   "+data.whiteListPII);
					initData(data);
					//alert('data '+json); 
				}				
			},
			 statusCode:function(status) {
			  console.log( "page not found"+ status );
			},
			error: function(e){
				if(e.status == 200){
					console.log("status:"+e.status+' statusText:'+e.statusText);
				}else if(e.status == 404){
					console.log("status:"+e.status+' statusText:'+e.statusText);
				}
			}
		});	
	}		
}); // on dom ready	


function getTimeStamp() {
	var now = new Date();
	return ((now.getMonth() + 1)+''+(now.getDate())+''+ now.getFullYear()+''+ now.getHours()+''+ ((now.getMinutes() < 10) ? ("0" + now.getMinutes()) : (now.getMinutes()))+''+ ((now.getSeconds() < 10) ? ("0" + now.getSeconds()) : (now.getSeconds())));
}

//Extarcting the Cytoscape network data from the ajax response
function initData(srcData){
	//var isEmptyNetwork = false;
	if(srcData.hasOwnProperty('cyjsFileName')){
		cyjsFileName = JSON.stringify(srcData.cyjsFileName);
	}
	
	if(srcData.hasOwnProperty('jsonFileName')){
		jsonFileName = JSON.stringify(srcData.jsonFileName);
	}
	
	if(srcData.hasOwnProperty('validCytoscapeNetwork')){
		validCytoscapeNetwork = JSON.stringify(srcData.validCytoscapeNetwork);
	}
	
	if(srcData.hasOwnProperty('cyjsFileContent')){
		cyjsFileContent = JSON.stringify(srcData.cyjsFileContent);
	}
	
	if(srcData.hasOwnProperty('jsonFileContent')){
		jsonFileContent = JSON.stringify(srcData.jsonFileContent);
	}
	
	var networkList = getNetworkList();
	var validCytoscapeNetworkList = getValidCytoscapeNetworkList();
	var cyDataList = getCyData();
	var cyStyleList = getCyStyle();
	
	/*if(cyDataList[0].elements.nodes.length == 0 && cyDataList[0].elements.edges.length == 0){
		alert('ZERO nodes and ZERO edges');
		isEmptyNetwork = true;
	}else{
		alert(cyDataList[0].elements.nodes.length + ' nodes and '+cyDataList[0].elements.edges.length + ' edges');
	}*/
	//if(!isEmptyNetwork){
		if(cyDataList.length>0 && cyStyleList.length>0){
			if(networkList && networkList.length> 0 && networkList.length > currentIndex){
				
				//Create network dropdown list				
				//populateNetworkList(cyDataList);
				populateNetworkList(networkList);
				
				//Add Zoom slider controller into the App
				createZoomSlider();				
				
				//Initialize all the Events handler
				initializeEvents();
				
				if(!isFullScreen){
					/*var obj = {};
					obj.event = 'cytoscapeLoadedEvent';
					obj.visible = true;
					var str = JSON.stringify(obj);
					parent.postMessage(str, "*");
					*/
					
					//This methode track the gadget impression 
					logSCDisplay();					
					
					dispatchPostMessageEvent('cytoscapeSuccessEvent');					
					//dispatchPostMessageEvent('cytoscapeLoadedEvent');	
				}
				
				isAppLoaded = true;
			}
		}
	//}
}

function getCyData(){
	return validateCyContent(cyjsFileContent);
}

function getCyStyle(){
	return validateCyContent(jsonFileContent);
}

function getNetworkList(){
	return validateCyContent(cyjsFileName);
}

function getValidCytoscapeNetworkList(){
	return validateCyContent(validCytoscapeNetwork);
}

function validateCyContent(objStr){
	var obj = JSON.parse(objStr);
	var isArr = Array.isArray(obj);
	if(isArr){
			return obj;
	}else{
		if(obj != null){
			var tempArr = [obj];
			return tempArr;
		}else{
			return [];
		}
	}
}

function createCyGraph(index){
	if(index >= 0){
		var cyDataArr = getCyData();
		var styleArr = getCyStyle();
		var validCytoscapeNetworkArr = getValidCytoscapeNetworkList();
		var cyStyleList = [];
		if(cyDataArr.length > 0 && cyDataArr.length > index){
			if(styleArr.length>0){
				if(Array.isArray(styleArr[index]) && styleArr[index].length > 0){
					var tempList = styleArr[index];					
					for(var i=0;i<tempList.length;i++){
						var tempObj = tempList[i];
						if(tempObj.hasOwnProperty('title') && tempObj.hasOwnProperty('style')){						
							if(tempObj.title == 'default'){							
								cyStyleList = tempObj.style;
								break;
							}
						}
					}
				}
			}	
			if(validCytoscapeNetworkArr[index]){
				//Cytoscape default data graph renderring
				//alert('valid Cytoscape Network');
				initCy( cyDataArr[index].elements, cyStyleList ); // cyStyleList[index][0].style
			}else{
				//alert('invalid cytoscape network');				
				$(".cyContainer").css('display','none');
				
				/*$(".preloader").css('display','none');
				$("#cy").css('display','none'); */
				$(".checkLabelDiv").css('display','none');
				$(".zoomDiv").css('display','none');
				$(".expIcon").css('display','none');
				$(".dwnImgDiv").css('display','none');				
				$("#dwnImgLink").css('display','none');
				$(".dwnNetDiv").css('display','inline-block');				
				$("#invalidNetwork").css('display','block');				
				updateGadgetHeight();
			}
		}
	}	
}	

function initCy( dataElements , dataStyle) {
	//Create the instance for Cytoscape Object 
	if(cyElement){
		cyElement.destroy();
	}
	
	var cyDiv = document.createElement('div');
	cyDiv.id = 'cy';
	$('.cyContainer').css('display','block');	
	$('.cyContainer').append(cyDiv);	
	
	cyElement = cytoscape({
		container: $('#cy')[0],
		elements: dataElements,		
		userZoomingEnabled: false,
		style: dataStyle,
		layout: { name: 'preset'}
	});
	var eles = cyElement.elements()
	if(eles.empty()){
		alert('Empty network');
	}else{
		//alert('non-empty network');
	}
	$('#cy').attr('tabIndex',1);	
	$(".preloader").css('display','none');
	$("#invalidNetwork").css('display','none');
	$(".dwnNetDiv").css('display','none');
	$(".checkLabelDiv").css('display','block');
	$(".zoomDiv").css('display','block');
	$(".expIcon").css('display','block');
	$(".dwnImgDiv").css('display','inline-block');
	$("#dwnImgLink").css('display','inline-block');
	
	
	cyElement.isEdgeChecked = true;
	cyElement.isNodeChecked = true;
	
	cyElement.minZoom(0);
	cyElement.maxZoom(10);
	
	cyElement.ready( function(){		
		if(zoomSlider){
			var t = cyElement.zoom() * 10;
			zoomSlider.setSliderValue(t); //cyElement.zoom();
		}		
	});	
	
	//This method will invoked on cytoscape node/edge selection
	cyElement.on('tap', 'node,edge', function(e){		
		var elementLabel = this.style().content;
		currentElement = this;
		isNodeEdge = true;		
		var isEdge = this.isEdge();
		var isNode = this.isNode();

		if(isNode){
			trackSCClick("Node label click");
		}else{
			trackSCClick("Edge label click");
		}
		
		if(elementLabel !='' && ((cyElement.isEdgeChecked && isEdge) ||  (cyElement.isNodeChecked && isNode))){
			var isQtip = this.hasOwnProperty('isQtipInjected');			
			if(!isQtip){				
				this['isQtipInjected'] = true;
				 //This method will associate tool-tip with selected nodes/edges.
				 this.qtip({						
						content: {
							text: function(api){
								return elementLabel;
							}
						},                                           
						position: {
							my: 'bottom center',
							at: 'top center'
						},
						style: {
							classes: 'qtip-bootstrap',
							tip: {
							width: 16,
							height: 8
							}
						}
					});
					//$(this).trigger(e);
					if(isEdge){
						var cy = this.cy();
						var container = cy.container();
						var cOff = container.getBoundingClientRect();
						var pos = e.cyRenderedPosition;						
						this._private.scratch.qtip.api.set('position.adjust.x', cOff.left + pos.x + window.pageXOffset);
						this._private.scratch.qtip.api.set('position.adjust.y', cOff.top + pos.y + window.pageYOffset);
					}
					currentElement._private.scratch.qtip.api.show();
			}
		}
		$('#cy').blur();
	});
	
	//This method will invoked on cytoscape graph click
	cyElement.on('tap', function(e){	
		//Remove focus from the ZOOM slider dragger elment
		$('#cy').focus();
		zoomSlider.removeFocus();
		if(isNodeEdge){				
			isNodeEdge = false;	
		}else{
			currentElement = null;
			var selectedElements = cyElement.$(':selected');
			
			if(selectedElements && selectedElements.length>0){
				for( var i = 0; i < selectedElements.length; i++ ){
					var node = selectedElements[i];
					node.unselect();
				}
			}
		}
	});
	
	//This method will invoked on mouse scroll on top of the Cytoscape graph area
	$('#cy').mousewheel(function(event, delta) {
		event.preventDefault();
		if (delta > 0){
			mouseScrollUp();  
		}else if (delta < 0){
			mouseScrollDown();
		}
	});
	
	//This method will do the cytoscape graph panning using keyboard down event
	$( '#cy' ).keydown(function( event ) {
		//event.stopPropagation();		
		var kc = event.keyCode;			
		var active = document.activeElement;
		var cyid = document.getElementById("cy");			
		if(active == cyid){
			var kc = event.keyCode;
			var obj = cyElement.pan();				
			var increment = 2;
			switch (kc) {
				case 38:// up
					obj.y = obj.y-increment;
				break;
				case 39:// right
					obj.x = obj.x+increment;						
				break;
				case 37:// left
					obj.x = obj.x-increment;
				break;
				case 40:// down
					obj.y = obj.y+increment;								
				break;
			}				
			cyElement.pan(obj);
		}
	});
	
	updateGadgetHeight();
}	

function initializeEvents(){
	//This method will invoked by reset button click event
	$('#resetLink').click(function(){
		if(currentIndex >= 0){
			trackSCClick("Reset click");
			createCyGraph(currentIndex,false);
			resetElements(true);
		}
	});
	
	//This method will invoked by node label checkbox change event
	$('#nodeLabel').change(function(){
		if(currentIndex >= 0){
			trackSCClick("Node label show/hide click");
			if(this.checked){
				cyElement.isNodeChecked = true;
				cyElement.$('node').css({'text-opacity': 1});
			}else{
				cyElement.isNodeChecked = false;
				cyElement.$('node').css({'text-opacity': 0});
			}	
		}				
	});
	
	//This method will invoked by edge label checkbox change event
	$("#edgeLabel").change(function(){
		if(currentIndex >= 0){
			trackSCClick("Edge label show/hide click");
			if(this.checked){
				cyElement.isEdgeChecked = true;
				cyElement.$('edge').css({'text-opacity': 1});
			}else{
				cyElement.isEdgeChecked = false;
				cyElement.$('edge').css({'text-opacity': 0});
			}
		}			
	});
	
	//This method will invoked by download-network label click event
	$("#dwnNetLink,#dwnNetId").click(function(){		
		if(currentIndex >= 0){
			trackSCClick("Download network file click");
			var indexStr = "&index="+ currentIndex;
			var cydownlaodURL;
			cydownlaodURL = serverURL[environment] + appName + "/CytoscapeDownloadZipServlet?doi=" + doi + indexStr;				
			//cydownlaodURL = "http://138.12.204.127:8080/CytoscapeWeb/CytoscapeDownloadZipServlet?doi=10.1016/S9999-9994%2815%2920790-8&index=1";			
			$(this).attr('href', cydownlaodURL);
		}			
	});	
	
	//This method will invoked by download-Link click event
	$("#dwnImgLink").click(function(){
		if(currentIndex >= 0){
			trackSCClick("Download image click");
			var ie = isIE ();		
			var options = {
				full: true,
				scale: 1
			};
			
			//The "download" attribute will support only HTML 5 modern browser except IE browser, so that we have handled seperate implementation using form submit for IE browser below IF condition.
			if(ie){
				var servletUrl = serverURL[environment] + appName +"/CytoscapeDownloadImageServlet?doi="+ doi;
				if(ie >= 11){
					//var url = "http://cert-elsevier-apps.sciverse.com/CytoscapeWeb/CytoscapeDownloadImageServlet?doi="+doi;
					$("#imgBase64Content").val(cyElement.png(options));
					$("#dummyForm").attr("action",servletUrl);
					$("#dummyForm").attr("method","POST");
					$("#dummyForm").submit();			
				}
			}else{
				if(cyElement){
					var baseimage = cyElement.png( options );
					var filename = $('#network_cmb option:selected').text()+".png";
					$(this).attr("download", filename);			
					$(this).attr('href', baseimage);
				}
			}			
		}
			
	});	
	
	//This method will invoked by download-Image click event
	$("#dwnImgId").click(function(){
		if(currentIndex >= 0){
			trackSCClick("Download image click");
			var ie = isIE ();
			var options = {
				full: true,
				scale: 1
			};
			
			if(ie){
				var servletUrl = serverURL[environment] + appName +"/CytoscapeDownloadImageServlet?doi="+ doi;
				if(ie >= 11){
					$("#imgBase64Content").val(cyElement.png(options));
					$("#dummyForm").attr("action",servletUrl);
					$("#dummyForm").attr("method","POST");
					$("#dummyForm").submit();			
				}
			}else{
				var baseimage = cyElement.png( options );
				var filename = $('#network_cmb option:selected').text()+".png";
				$(this).attr("download", filename);
				$(this).attr('href', baseimage);	
			}
		}		
	});	
	
	//This method will invoked by newwideow/Explore-link click event
	$("#newWindowIcon, #explrImg, #explrLink").click(function(evt){	
		openNewWindow(evt.target);
	});
	
	$('.cytoscapeImg').click(function(){
		//trackSCClick('Cytoscape Logo clicks');
		window.open('http://www.cytoscape.org/','_blank');
	});
}

//This method will open the Article content page with new window 
function openArticlePage(){
	trackSCClick("More information click");
	var url = 'http://www.elsevier.com/about/content-innovation/cytoscape';		
	window.open(url,'_blank');
}

//This method will open the full screen view with new window 
function openNewWindow(obj){	
	if(currentIndex >=0){
		if(obj.id === "newWindowIcon"){
			trackSCClick("Detachable icon click");
		}else if(obj.id === "explrImg"){
			trackSCClick("Explore icon click");
		}else{
			trackSCClick("Explore link click");
		}
		//isFullScreen = true;
		var myObject = {
			index: currentIndex,
			doi: doi,
			domain: environment,
			isFullScreen: true
		};
		var url;
		url = serverURL[environment]+ appName+"/index.html?"+$.param( myObject );		
		window.open(url,'_blank');	
	}
}

//This method will help to validated the IE browser version detect 
function isIE() {
	var rv = -1;
	if (navigator.appName == 'Microsoft Internet Explorer'){
		var ua = navigator.userAgent;
		var re  = new RegExp("MSIE ([0-9]{1,}[\.0-9]{0,})");
		if (re.exec(ua) != null)
		rv = parseFloat( RegExp.$1 );
	}else if (navigator.appName == 'Netscape'){
		var ua = navigator.userAgent;
		var re  = new RegExp("Trident/.*rv:([0-9]{1,}[\.0-9]{0,})");
		if (re.exec(ua) != null){
			rv = parseFloat( RegExp.$1 );
		}
	}
	return (rv != -1)?rv:false;
}		

//This method will invode on reset elements click event
function resetElements(isNetworkChange){		
	$('#nodeLabel').prop('checked',true);
	$('#edgeLabel').prop('checked',true);	
}

function mouseScrollUp(){	
	if(zoomSlider){
		var t = cyElement.zoom() * 10;
		(t<100)?t++:t=100;
		zoomSlider.setSliderValue(t);
	}	
}

function mouseScrollDown(){	
	if(zoomSlider){
		var t = cyElement.zoom() * 10;
		(t>0)?t--:t=0;
		zoomSlider.setSliderValue(t);
	}
}

//This method will create the "rangeSlider" instance and append into the app DOM elements
function createZoomSlider(){
	// Horizontal range slider
	zoomSlider = new rangeSlider(document.getElementById('range-Slider'), {		
		drag: function(v) {           
			networkGraphZoom(v);
		}
	});
	
	$( ".range-slider" ).on( "mousedown", function() {
		trackSCClick("Zoom slider click");
	});
}

//This method will invoked by range slider thumb dragding event
function networkGraphZoom(val){
	if(cyElement){
		var tval = val/10;
		cyElement.zoom(tval);			
		var selectedElements = cyElement.$(':selected');		
		if(selectedElements && selectedElements.length>0){		
			var node = currentElement;
			cyElement.center(node);	
		}else{
			cyElement.center();	
		}	
	}		
}

//This method will populate network lists on the dropdown selection
function populateNetworkList(networkList){	
	if(Array.isArray(networkList) && networkList.length > 0 && networkList.length > currentIndex){
		for(var i = 0; i < networkList.length; i++) {		
			var opt = networkList[i];	
			var network = JSON.stringify(opt);
			var el = document.createElement("option");
			/*if(opt.hasOwnProperty('data')){
				el.textContent = opt.data.name;
			}*/
			if(network.lastIndexOf("/")!=-1){
				el.textContent = network.slice(network.lastIndexOf("/")+1,network.lastIndexOf(".cyjs"));
			}else{
				el.textContent = network.slice(1,network.lastIndexOf(".cyjs"));
			}
			el.value = i;
			$('#network_cmb').append(el);
		}
		
		$('#network_cmb').change(function(){
			trackSCClick("Network selection click");
			currentIndex = $(this).val();			
			createCyGraph(currentIndex);
			resetElements(false);
		});
		
		if(isFullScreen){	
			if(currentIndex){
				$('#network_cmb').val(currentIndex);
			}
		}
		$("#network_cmb").change();			
	}
}


function getQueryVariable(){
	var queryString;
	if (window.location.search.split('?').length > 1) {
		var params = window.location.search.split('?')[1].split('&');
		if(params.length>0){
			queryString= {};
			for (var i = 0; i < params.length; i++) {
				var key = params[i].split('=')[0];
				var value = params[i].split('=')[1];
				var dc = decodeURIComponent(value);				
				queryString[key] = decodeURIComponent(value);
			}
		}
	}			
	return (queryString)?queryString:false;		
}

//STEP #6 - To track all the Site-Catalyst SC click events
function trackSCClick(evtLabel){	
	if(isFullScreen){
		evtLabel = "Full View - "+evtLabel
	}else{
		evtLabel = "Gadget View - "+evtLabel
	}
	scGadgetClickEvent = evtLabel;
	logSCUsage();
}

function updateGadgetHeight(){
	dispatchPostMessageEvent('cytoscapeHeightChangeEvent');	
}


//Dispath postMessage event to the parent html page 
function dispatchPostMessageEvent(event , isVisible){
	var obj = {};
	obj.eventName = event;
	obj.isVisible = (typeof isVisible !== 'undefined') ? isVisible : true;
	obj.appHeight = $('#appContainerId').outerHeight(true);
	var str = JSON.stringify(obj);
	parent.postMessage(str, "*");	
}

/***   Cytoscape code END  ***/