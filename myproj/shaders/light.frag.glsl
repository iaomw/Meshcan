#version 330 core

in vec4 myvertex;
in vec3 mynormal;

out vec4 finalcolor;

uniform mat4 myview_matrix;
uniform mat3 mynormal_matrix;

uniform int type;

uniform vec4 kd;

vec4 mylight_position = vec4(0,0,4,1);
vec4 mylight_color = vec4(1,1,1,0);

void main (void) 
{        
	vec3 eyepos = vec3(0,0,0) ; 

	vec4 _mypos = myview_matrix * myvertex ; 
    vec3 mypos = _mypos.xyz / _mypos.w ;

	vec4 _lightpos = mylight_position;
   	vec3 lightpos = _lightpos.xyz / _lightpos.w;

	vec3 normal = normalize(mynormal_matrix*mynormal);

	vec3 eyedir = normalize(eyepos - mypos) ;
	vec3 lightdir = normalize (lightpos - mypos) ;

	if (dot(normal, eyedir) < -0.1) normal = -normal;

	vec3 reflectdir = normalize( reflect(-lightdir, normal) );

    finalcolor =  vec4(0.1,0.1,0.1,0) + 
					kd * mylight_color * max( dot(lightdir, normal), 0.0) +
					kd * mylight_color * pow( max( dot(reflectdir, eyedir), 0.0 ), 60 );
    if (type == 1) finalcolor = kd;
}