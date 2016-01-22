/*
 *
 * =======================================================================
 *  THE JAVASCRIPT CUSTOM SLIDER CONTROLLER
 *  Author: GOPIKANNAN	
 * =======================================================================
 *
 */

/*
* Sample code for Horizontal and Vertical range slider instances adding on the HTML JS code
*	// Horizontal range slider
*	var hSlider = new rangeSlider(document.getElementById('range-slider-1'), {
		 width: 200,
		 height: 20,
*        drag: function(v) {
*            testResult.children[0].innerHTML = v + '%';
*        }
*    });
*
*	
*	// Vertical range slider
*	var vSlider = new rangeSlider(document.getElementById('range-slider-2'), {
*       value: 20,
*       vertical: true,
*       drag: function(v) {
*           testResult.children[1].innerHTML = v + '%';	
*			temp(v);			
*       },		
*   });
*
*/
function rangeSlider(elem, config) {

    var html = document.documentElement,
        range = document.createElement('div'),
        dragger = document.createElement('span'),
        down = false,
        rangeWidth, rangeOffset, draggerWidth, cachePosition;
		
	var minWidth = 50;
	var minHeight = 50;

    var defaults = {
        value: 0, // set default value on initiation from `0` to `100` (percentage based)
        vertical: false, // vertical or horizontal?
        rangeClass: "", // add extra custom class for the range slider track
        draggerClass: "", // add extra custom class for the range slider dragger
		width:"",
		height:"",
        drag: function(v) { /* console.log(v); */ } // function to return the range slider value into something		
    };

    for (var i in defaults) {
        if (typeof config[i] == "undefined") config[i] = defaults[i];
    }

    function addEventTo(el, ev, fn) {
        if (el.addEventListener) {
            el.addEventListener(ev, fn, false);
        } else if (el.attachEvent) {
            el.attachEvent('on' + ev, fn);
        } else {
            el['on' + ev] = fn;
        }
    }

    var isVertical = config.vertical;

    elem.className = (elem.className + ' range-slider ' + (isVertical ? 'range-slider-vertical' : 'range-slider-horizontal')).replace(/^ +/, "");
	if(config.width!="" && !isNaN(config.width)){
		var cols = document.getElementsByClassName('range-slider');
		
		for(i=0; i<cols.length; i++) {
			cols[i].style.width = ((Number(config.width) > minWidth)?config.width:minWidth) + "px";
		  }	
	}		
	if(config.height!="" && !isNaN(config.height) && config.height > minHeight ){
		var cols = document.getElementsByClassName('range-slider');
		for(i=0; i<cols.length; i++) {
			cols[i].style.height = config.height + "px";
		  }	
	}
	
    range.className = ('range-slider-track ' + config.rangeClass).replace(/ +$/, "");
    dragger.className = ('dragger ' + config.draggerClass).replace(/ +$/, "");	
	dragger.tabIndex  = 0;
	range.tabIndex  = 0;

    addEventTo(range, "mousedown", function(e) {
        html.className = (html.className + ' no-select').replace(/^ +/, "");
		updateSliderElements();
        down = true;
		updateDragger(e);				
        return false;
    });

    addEventTo(document, "mousemove", function(e) {
        updateDragger(e);
    });

    addEventTo(document, "mouseup", function(e) {
        html.className = html.className.replace(/(^| )no-select( |$)/g, "");
        down = false;		
    });

    addEventTo(window, "resize", function(e) {		
		dragger.blur();
        var woh = dragger[!isVertical ? 'offsetWidth' : 'offsetHeight'];
        dragger.style[!isVertical ? 'left' : 'top'] = (((cachePosition / 100) * range[!isVertical ? 'offsetWidth' : 'offsetHeight']) - (woh / 2)) + 'px';
		down = false;
    });
	
	addEventTo(document.documentElement, "keydown", function(e) {		
		var active = document.activeElement;		
		if(active == dragger){
			e.preventDefault();			
			var pVal = getPercentage() ;//this.getSliderValue;
			//console.log(pVal);
			var isKeyAllow = false;			
			var kc = e.keyCode;				
			switch (kc) {
				case 38:	// up
				case 39:	// right
					isKeyAllow = true;
					pVal++										
					break;
				case 37:	// left
				case 40:	// down
					isKeyAllow = true;
					pVal--								
					break;
			}
			
			if(isKeyAllow){								
				if(pVal <= 0){
					pVal = 0;
				}else if(pVal >= 100){
					pVal = 100;
				}				
				dragger.style[!isVertical ? 'left' : 'top'] = (getXPosition(pVal)) + 'px';
				cachePosition = pVal;
				config.drag(pVal);											
			}
		}else{			
			//console.log(e.keyCode+" drag not focus ");
		}		
    });
	
	addEventTo(range, "focus", function(e) {		
		dragger.focus();
	});
	
	addEventTo(range, "blur", function(e) {
		dragger.blur();
	});
	
		
	function updateSliderElements(){
		rangeWidth = range[!isVertical ? 'offsetWidth' : 'offsetHeight'];
        rangeOffset = range[!isVertical ? 'offsetLeft' : 'offsetTop'];
        draggerWidth = dragger[!isVertical ? 'offsetWidth' : 'offsetHeight'];		
	}
	
    function updateDragger(e) {
        e = e || window.event;
        var pos = !isVertical ? e.pageX : e.pageY;
        if (!pos) {
            pos = !isVertical ? e.clientX + document.body.scrollLeft + document.documentElement.scrollLeft : e.clientY + document.body.scrollTop + document.documentElement.scrollTop;
        }
        if (down && pos >= rangeOffset && pos <= (rangeOffset + rangeWidth)) {
            dragger.style[!isVertical ? 'left' : 'top'] = (pos - rangeOffset - (draggerWidth / 2)) + 'px';
            cachePosition = Math.round(((pos - rangeOffset) / rangeWidth) * 100);			
            config.drag(cachePosition);
        }
    }
	
	this.removeFocus = function(){
		dragger.blur();
	}
	
	this.setFocus = function(){
		dragger.focus();
	}
	
	this.setSliderValue = function(val){
		cachePosition = val;
		dragger.style[!isVertical ? 'left' : 'top'] = (getXPosition(val)) + 'px';		
		config.drag(val);
	}
	
	this.getSliderValue = function(){	
		return getPercentage();
	}
	
	function getPercentage(){		
		var woh = dragger[!isVertical ? 'offsetWidth' : 'offsetHeight'];			
		var dX = dragger[!isVertical ? 'offsetLeft' : 'offsetTop'] + (woh / 2);
		var v = Math.round((dX/rangeWidth)*100);
		return v;
	}
	
	function getXPosition(val){
		var woh = dragger[!isVertical ? 'offsetWidth' : 'offsetHeight'];
		var xpos = ((val / 100) * range[!isVertical ? 'offsetWidth' : 'offsetHeight']) - (woh / 2);
		return xpos;
	}

    function initDragger() {
        var woh = dragger[!isVertical ? 'offsetWidth' : 'offsetHeight'];
        cachePosition = ((config.value / 100) * range[!isVertical ? 'offsetWidth' : 'offsetHeight']);
        dragger.style[!isVertical ? 'left' : 'top'] = (cachePosition - (woh / 2)) + 'px';
        config.drag(config.value);		
    }

    range.appendChild(dragger);
    elem.appendChild(range);
    initDragger();	
}