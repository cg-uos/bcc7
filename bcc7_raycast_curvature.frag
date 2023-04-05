#version 330 core

in vec2 vTexCoord;

uniform sampler2D	tex_back;
uniform sampler2D	tex_front;
uniform sampler3D	tex_volume;
uniform sampler2D   tex_colormap;
uniform sampler2D   tex_colormap_2d;
uniform vec3        scale_norm;
uniform vec3        offset_norm;
uniform vec3        scale_lattice;
uniform vec3        offset_lattice;
uniform	vec3		dim;
uniform float       dim_max;
uniform float		level;
uniform mat4        MV;
uniform float       scale_step;
uniform float		scale_delta;

out vec4 fColor;

#define	DISCARD					\
	discard;




//#define SCALE_M     (0.5*0.00260416666)
//#define SCALE_Tj    (0.5*0.02083333333)
//#define SCALE_Tij   (0.5*0.125)
#define SCALE_M     (1.0/192.0)
#define SCALE_Tj    (1.0/24.0)
#define SCALE_Tij   (1.0/4.0)

float	c[30];
vec3    xtet;
vec4    u;
ivec3	org;

int[4]    vecR;
int[8]    vecP;
ivec3	type_R;

int	idx;
int	idx_R, idx_P;


#define R000    vecR[0]
#define R011    vecR[1]
#define R101    vecR[2]
#define R110    vecR[3]

#define P000    vecP[0]
#define P001    vecP[1]
#define P100    vecP[2]
#define P110    vecP[3]
#define P111    vecP[4]
#define P011    vecP[5]


#if	1
#define	GET_DATA(texcoords)	texelFetch(tex_volume, texcoords, 0).r
#else
int offset_x;
#define	GET_DATA(texcoords)	texelFetch(tex_volume, 	\
		ivec3(	(texcoords.x>>1) + offset_x*(texcoords.x&1),	\
				(texcoords.y>>1),					\
				(texcoords.z>>1)), 0).r
#endif



#define	c0	c[0]
#define	c1	c[1]
#define	c2	c[2]
#define	c3	c[3]
#define	c4	c[4]
#define	c5	c[5]
#define	c6	c[6]
#define	c7	c[7]
#define	c8	c[8]
#define	c9	c[9]
#define	c10	c[10]
#define	c11	c[11]
#define	c12	c[12]
#define	c13	c[13]
#define	c14	c[14]
#define	c15	c[15]
#define	c16	c[16]
#define	c17	c[17]
#define	c18	c[18]
#define	c19	c[19]
#define	c20	c[20]
#define	c21	c[21]
#define	c22	c[22]
#define	c23	c[23]
#define	c24	c[24]
#define	c25	c[25]
#define	c26	c[26]
#define	c27	c[27]
#define	c28	c[28]
#define	c29	c[29]


void preprocess(vec3 p_in)
{
	ivec3	lower = ivec3(p_in);

	type_R = ivec3(lower.x&1, lower.y&1, lower.z&1);

//	int	parity = type_R.x ^ type_R.y ^ type_R.z;
	int	parity = (type_R.x + type_R.y + type_R.z)&1;

//	type_R = parity*(type_R ^ ivec3(1)) + (1 ^ parity)*type_R;
	type_R = parity*(ivec3(1,1,1) - type_R) + (1-parity)*type_R;

    org = lower + type_R;


	vec3	x_cube = (p_in - vec3(org)) * vec3(1 - 2*type_R);

    ivec3    type_P = ivec3(x_cube.y >= x_cube.x, x_cube.z >= x_cube.x, x_cube.z >= x_cube.y);


    idx_R = (type_R.x<<1) + type_R.y;
    idx_P = (1-type_P.x)*(1-type_P.y)*type_P.z  // (0,0,0) or (0,0,1)
            + type_P.x*(1-type_P.z)*(2+type_P.y)    // (1,0,0) or (1,1,0)
            + type_P.y*type_P.z*(5-type_P.x);   // (1,1,1) or (0,1,1)

	idx = idx_R*6 + idx_P;

    vecR = int[](int(idx_R==0),int(idx_R==1),int(idx_R==2),int(idx_R==3));
    vecP = int[](int(idx_P==0),int(idx_P==1),int(idx_P==2),int(idx_P==3),int(idx_P==4),int(idx_P==5),-1,-1);


    xtet = float(P000)*x_cube.xyz
        +  float(P001)*x_cube.xzy
        +  float(P100)*x_cube.yxz
        +  float(P110)*x_cube.yzx
        +  float(P111)*x_cube.zyx
        +  float(P011)*x_cube.zxy;

	u = vec4(1-xtet.x, xtet.x-xtet.y, xtet.y-xtet.z, xtet.z);
}

void fetch_coefficients(void)
{
	ivec3	o = org;
    ivec3	sgn = 2-4*ivec3(R101+R110, R011+R110, R011+R101);
    ivec3	R[3] = ivec3[](ivec3((P000+P001)*sgn.x, (P100+P110)*sgn.y, (P111+P011)*sgn.z),
                            ivec3((P100+P011)*sgn.x, (P000+P111)*sgn.y, (P001+P110)*sgn.z),
                            ivec3((P110+P111)*sgn.x, (P001+P011)*sgn.y, (P000+P100)*sgn.z));


#define	FETCH(var)	var = GET_DATA(o);


					FETCH(c0); o -= R[0];		FETCH(c10); o += (R[0]<<1);	FETCH(c9);
	o -= R[1];		FETCH(c20); o -= R[0];		FETCH(c12); o += (R[1]<<1);	FETCH(c11);
	o += R[0];		FETCH(c19); o -= R[2];		FETCH(c29); o -= R[1];		FETCH(c18);
	o -= R[0];		FETCH(c14); o += R[1];		FETCH(c16); o += (R[2]<<1);	FETCH(c15);
	o += R[0];		FETCH(c28); o -= R[1];		FETCH(c17); o -= R[0];		FETCH(c13);
	o += ((-R[0]+R[1]-3*R[2])>>1);	FETCH(c6);
	o -= R[1];		FETCH(c8); o += R[0];		FETCH(c4); o += (R[1]<<1);	FETCH(c26);
	o -= R[1];		FETCH(c2); o += (R[2]<<1);	FETCH(c27); o -= R[2];		FETCH(c1);
	o += R[1];		FETCH(c25); o -= (R[1]<<1);	FETCH(c3); o -= R[0];		FETCH(c7);
	o += R[1];		FETCH(c5); o += (R[0]<<1);	FETCH(c21); o -= R[2];		FETCH(c22);
	o -= R[1];		FETCH(c24); o += R[2];		FETCH(c23); 
#undef  FETCH
}
#define	EVAL(p)	(preprocess(p), fetch_coefficients(), eval_M())
float eval_M(void)
{
#define	x	(xtet.r)
#define	y	(xtet.g)
#define	z	(xtet.b)

#define	SQR(X)	((X)*(X))
#define	TRI(X)	((X)*(X)*(X))
#define	QUAD(X)	((X)*(X)*(X)*(X))

#define	xpy	(x+y)
#define	xmy	(x-y)
#define	xpz	(x+z)
#define	xmz	(x-z)
#define	ypz	(y+z)
#define	ymz	(y-z)

	float	v8 = 0.0;
	float	v4 = 0.0;
	float	v2 = 0.0;
	float	t1;
	t1 = SQR(x);
	v8 -= SQR(t1)*(c9+c19);
	t1 = SQR(x-1.0);
	v8 += SQR(t1)*(c0 - c5 - c6 - c7 - c8 + c10 + c11 + c12 + c14);
	t1 = SQR(x-3.0);
	v8 += SQR(t1)*c0;
	t1 = SQR(x-2.0);
	v8 -= SQR(t1)*(c0 - c5 + c11);
	t1 = SQR(x+1.0);
	v8 += SQR(t1)*c9;
	t1 = SQR(z);
	v8 += SQR(t1)*(c1 + c3 - c13 - c15 - c17 + c21 + c25 + c27 - c28);
	t1 = SQR(z-1.0);
	v8 -= SQR(t1)* (c2 + c4);
	t1 = SQR(z+1.0);
	v8 -= SQR(t1)*(c1 + c3 - c13);
	t1 = SQR(z-2.0);
	v8 += SQR(t1)*c2;
	t1 = SQR(z+2.0);
	v8 += SQR(t1)*c1;

	v4 -= TRI(xpy)*(2.0+xmy)*(c9 - c19);
	t1 = xmy-2.0;
	v4 -= TRI(t1)*(xpy-4.0)*(c0-c11);
	t1 = ypz-2.0;
	v4 += TRI(t1)*(ymz+2.0)*(c2-c4);
	t1 = ymz-2.0;
	v4 += TRI(t1)*(ypz+2.0)*(c1-c3);

	v2 -= TRI(xmy)*(
			2.0*(
				2.0*(c9 - c20)
				- z*(c23 - c24)
				)
			- xpy*(2.0*(c19 - c20) + c23 + c24));
	v2 -= TRI(xpz)* (
			2.0*(
				2.0*(c1 + c3 + c9 + c13 - c17) +
   				y*(c1 - c3 - c9 + c19 - c21 )
				)
			- xmz*(2.0* (c13 - c17) + c1 + c3 - c9 - c19 + c21) );
	v2 += TRI(xmz)*(
			2.0*(
				-2.0*( c2 + c4 + c9 - c18)
				-y*(c2 - c4 - c9 + c19 - c22 )
				)
			+ xpz* (c2 + c4 - c9 - 2.0* c18 - c19 + c22) );
	v2 -= TRI(ypz)*(
			2.0*(
				2.0* (c1 - c3 + c13 - c15)
				- x*(c13 - c15 - c17 - c21 + c28)
				)
			+ ymz* (2.0* (c3 - c25) - c13 + c15 - c17 - c21 + c28) ) ;
	v2 -= TRI(ymz)*(
			2.0*(
				2.0* (c2 - c4 - c16)
				+ x*(c16 + c18 + c22 - c29)
				)
			- ypz* (-2.0* (c4 - c26) + c22 - c16 + c18 - c29));
	t1 = xpy-2.0;
	v2 += TRI(t1)*(
			2.0*(
				2.0* (c0 - c11)
				+ (c5 + c6 - c7 - c8 + c14)
				+ z* (c5 - c6 - c7 + c8 - c14 )
				)
			- xmy* (c5 + c6 - c7 - c8 + c14 - 2.0* ( c11 - c12)));
	t1 = xpz-2.0;
	v2 += TRI(t1)*(
			2.0*(
				3.0* c0 +
				(c2 + c4 + c11 - c14)
				+ 2.0* (c5 - c6)
				-  y* (c0 - c2 + c4 - c11 - c14)
				)
			- xmz* (c0 - c2 - c4 + c11 + c14 + 2.0* (c5 - c6)) );
	t1 = xmz-2.0;
	v2 += TRI(t1)* (
			2.0*(
				3.0* c0
				+  (c1 + c3 + c11) -
  				 y *(c0 - c1 + c3 - c11)
				)
			-xpz* (c0 - c1 - c3 + c11));

	return (4.0*v8 + 2.0*v4 + v2)*SCALE_M;

#undef	x
#undef	y
#undef	z

#undef	xpy
#undef	xmy
#undef	xpz
#undef	xmz
#undef	ypz
#undef	ymz


}


#define	u0	u.x
#define	u1	u.y
#define	u2	u.z
#define	u3	u.w

vec3 eval_Tj(void)
{
	float	u00 = u0*u0;
	float	u11 = u1*u1;
	float	u22 = u2*u2;
	float	u33 = u3*u3;

    float   u000 = u00*u0;

    float   u001 = u00*u1;
    float   u002 = u00*u2;
    float   u003 = u00*u3;

    float   u011 = u0*u11;
    float   u012 = u0*u1*u2;
    float   u013 = u0*u1*u3;
    float   u022 = u0*u22;
    float   u023 = u0*u2*u3;
    float   u033 = u0*u33;

    float   u111 = u11*u1;
    float   u112 = u11*u2;
    float   u113 = u11*u3;
    float   u122 = u1*u22;
    float   u123 = u1*u2*u3;
    float   u133 = u1*u33;

    float   u222 = u22*u2;
    float   u223 = u22*u3;
    float   u233 = u2*u33;
    
    float   u333 = u33*u3;

    vec3    expr;
	expr.x = 
		  u012*12*(6*(c11-c0) + 7*(c1+c2-c3-c4) + 2*(c19-c9) + c5 + c6 - c7 - c8 )
		+ u013*24*(3*(c11-c0) + 5*(c1-c3) + c19 + 2*(c2-c4) + c5 - c7 - c9) 
		+ u011*6*(7*(c1+c2-c3-c4) + 3*(c11 - c12) + c19 - c20 + c5 + c6 - c7 - c8) 
		+ u023*24*(5*(c11-c0) + 3*(c1-c3) - c13 + c15 + 2*(c19-c9) + c2 - c4 ) 
		+ u022*6*(10*(c11-c0) + 4*(c1+c19+c2-c3-c4-c9) - c13 - c14 + c15 + c16 )
		+ u033*24*(2*(c11-c0+c1-c3) - c13 + c15 + c19 - c9)
		+ u000*8*(c1 + c11 - c12 + c2 - c3 - c4 + c5 + c6 - c7 - c8) 
		+ u001*12*(3*(c1+c2-c3-c4) + 2*(c11 - c12) + c5 + c6 - c7 - c8) 
		+ u002*12*(4*(-c0+c11) + 3*(c1+c2-c3-c4) + c5 + c6 - c7 - c8) 
		+ u003*24*(2*(-c0 + c1 + c11-c3) + c2 - c4 + c5 - c7)
		+ u123*12*(7*(-c0+c11+c19-c9) + 6*(c1-c3) - c13 + c15 - c17 + 2*(c2-c4) + c28 )
		+ u122*3*(14*(-c0 +c11+c19-c9) + 8*(c1+c2-c3-c4) - c13 - c14 + c15 + c16 - c17 - c18 + c28 + c29 )
		+ u133*12*(3*(-c0+c11+c19-c9) + 4*(c1-c3) - c13 + c15 - c17 + c28 )
		+ u112*3*(8*(-c0+c11+c19-c9) + 14*(c1+c2-c3-c4) + c21 + c22 - c23 - c24 + c5 + c6 - c7 - c8 )
		+ u113*6*(4*(-c0+c11+c19+c2-c4-c9) + 10*(c1-c3) + c21 - c23 + c5 - c7 )
		+ u111*(14*(c1+c2-c3-c4) + 4*(c11 - c12 + c19 - c20) + c21 + c22 - c23 - c24 + c5 + c6 - c7 - c8) 
		+ u233*12*(3*(-c0 + c11 + c19 - c9) - c13 + c15 - c17 + 2*(c25-c3) + c28 )
		+ u223*6*(7*(-c0 + c11 + c19 - c9) - c13 + c15 - c17 + 3*(c25-c3) + c26 + c28 - c4 )
		+ u222*(14*(-c0 + c11+c19-c9) - c13 - c14 + c15 + c16 - c17 - c18 + 4*(c25 + c26 -c3-c4) + c28 + c29 )
		+ u333*8*(-c0 + c11 - c13 + c15 - c17 + c19 + c25 + c28 - c3 - c9);
	expr.y = 
		  u012*12*(10*(c1-c2) + 3*(c13 - c14 +c3-c4) + c17 - c18 + 2*(c5 - c6) ) 
		+ u013*24*(3*(-c0+c13) + 5*(c1-c2) + c17 + 2*(c3 - c4) + c5 - c6 - c9) 
		+ u011*6*(7*(c1-c2+c3-c4) + 3*(c13 - c14) + c17 - c18 + c5 - c6 + c7 - c8) 
		+ u023*24*(2*(-c0+c13) + 6*(c1-c2) - c11 + c15 + c17 + c3 - c4 + c5 - c6 - c9) 
		+ u022*6*(12*(c1-c2) + 2*(c13 - c14 + c3 - c4 + c5 - c6) + c15 - c16 + c17 - c18 )
		+ u033*24*(2*(-c0 + c1 +c13 - c2) - c11 + c15 + c17 - c9)
		+ u000*8*(c1 + c13 - c14 - c2 + c3 - c4 + c5 - c6 + c7 - c8) 
		+ u001*12*(3*(c1-c2+c3-c4) + 2*(c13 - c14) + c5 - c6 + c7 - c8) 
		+ u002*24*(2*(c1-c2) + c13 - c14  + c3 - c4 + c5 - c6) 
		+ u003*24*(2*(-c0 + c1 + c13 - c2) + c3 - c4 + c5 - c6)
		+ u123*12*(3*(-c0 +c13+c17-c9) + 12*(c1-c2) - c11 + c15 - c19 + c21 - c22 + c28 + 2*(c3 - c4) + c5 - c6 )
		+ u122*3*(24*(c1-c2) + 4*(c3 - c4) + 3*(c13 - c14 + c17 - c18 ) + 2*(c21 - c22 + c5 - c6) + c15 - c16 + c28 - c29 )
		+ u133*12*(3*(-c0+c13+c17-c9) + 4*(c1-c2) - c11 + c15 - c19 + c28 )
		+ u112*6*(10*(c1-c2) + 2*(c13 - c14 + c17 - c18 ) + c21 - c22 + 4*(c3 - c4) + c5 - c6) 
		+ u113*6*( 10*(c1 - c2) +4*(-c0 + c13 + c17 + c3 - c4 - c9) + c21 - c22 + c5 - c6 )
		+ u111*( 14*(c1 - c2 + c3 - c4) + 4*(c13 - c14 + c17 - c18) + c21 - c22 + c23 - c24 + c5 - c6 + c7 - c8) 
		+ u233*24*(-c0 + 2*(c1-c2) - c11 + c13 + c15 + c17 - c19 + c28 - c9) 
		+ u223*6*(  12*(c1 - c2) + 2*(-c0 - c11 + c13 + c15 + c17 - c19 + c28 - c9) + c21 - c22 + c25 - c26 + c3 - c4 + c5 - c6 )
		+ u222*2*(12*(c1-c2) + c13 - c14 + c15 - c16 + c17 - c18  + c21 - c22 + c25 - c26 + c28 - c29 + c3 - c4 + c5 - c6) 
		+ u333*8*(-c0 - c11 + c13 + c15 + c17 - c19 - c2 + c27 + c28 - c9);
	expr.z = 
		  u012*12*(12*(c9-c0) + 3*(c1+c2-c5-c6) + 2*(c19-c11) - c13 - c14 + c17 + c18 + c3 + c4 - c7 - c8 )
		+ u013*24*(6*(c9-c0) + 2*(c1-c5) - c11 - c13 + c17 + c19 + c2 + c3 - c6 - c7 )
		+ u011*6*(12*(c9-c0) + 2*(c1 + c2 + c3 + c4 - c5 - c6 - c7 - c8 ) - c11 - c12 - c13 - c14 + c17 + c18 + c19 + c20 )
		+ u023*24*(5*(c9-c0) + 3*(c1-c5) + 2*(c19-c11) - c13 + c17 + c2 - c6 )
		+ u022*6*(10*(c9-c0) + 4*(c1 - c11 + c19 + c2 - c5 - c6 ) - c13 - c14 + c17 + c18 )
		+ u033*24*(2*(-c0 + c1-c5+c9) - c11 - c13 + c17 + c19 )
		+ u000*8*(c1 - c10 + c2 + c3 + c4 - c5 - c6 - c7 - c8 + c9) 
		+ u001*24*(2*(c9-c0) + c1 + c2 + c3 + c4 - c5 - c6 - c7 - c8  ) 
		+ u002*12*(4*(c9-c0) + 3*(c1 + c2-c5-c6) + c3 + c4 - c7 - c8 )
		+ u003*24*(2*(-c0 + c1-c5+c9) + c2 + c3 - c6 - c7 )
		+ u123*12*(10*(c9-c0) + 4*(-c11+c19) + 2*(-c13 + c17 ) + 3*(c21-c5) + c22 - c6 )
		+ u122*6*(10*(c9-c0) + 4*(-c11+c19) - c13 - c14 + c17 + c18 + 2*(c21 + c22 - c5 - c6) )
		+ u133*24*(2*(c9-c0) - c11 - c13 + c17 + c19 + c21 - c5)
		+ u112*3*(24*(c9-c0) + 4*(c19-c11 ) + 3*(-c13 - c14 + c17 + c18 + c21 + c22 - c5 - c6 ) + c23 + c24 - c7 - c8 )
		+ u113*6*(12*(c9-c0) + 2*(-c11 - c13 + c17 + c19 + c21 - c5) + c22 + c23 - c6 - c7 )
		+ u111*2*(12*(c9-c0) - c11 - c12 - c13 - c14 + c17 + c18 + c19 + c20 + c21 + c22 + c23 + c24 - c5 - c6 - c7 - c8 )
		+ u233*12*(3*(c9-c0 - c11+c19) - c13 - c15 + c17  + 2*(c21-c5) + c28 )
		+ u223*6*(7*(c9-c0 - c11+c19) - c13 - c15 + c17 + 3*(c21-c5) + c22 + c28 - c6 )
		+ u222*(14*(c9-c0 - c11+c19) - c13 - c14 - c15 - c16 + c17 + c18 + 4*(c21 + c22-c5-c6) + c28 + c29 )
		+ u333*8*(-c0 - c11 - c13 - c15 + c17 + c19 + c21 + c28 - c5 + c9);	

    return expr*SCALE_Tj;
}

#undef	u0
#undef	u1
#undef	u2
#undef	u3



vec3 compute_gradient(vec3 p)
{
    vec3   expr = eval_Tj();
    vec3 g = float(P000)*expr.zxy
        +  float(P001)*expr.zyx
        +  float(P100)*expr.xzy
        +  float(P110)*expr.yzx
        +  float(P111)*expr.yxz
        +  float(P011)*expr.xyz;
    return 0.5*g*vec3(1-2*type_R)*scale_norm;
}

float[6] eval_Tij(void)
{
    float   expr[6];
    float   xx = u.x*u.x;
    float   yy = u.y*u.y;
    float   zz = u.z*u.z;
    float   ww = u.w*u.w;
    float   xy = u.x*u.y;
    float   xz = u.x*u.z;
    float   xw = u.x*u.w;
    float   yz = u.y*u.z;
    float   yw = u.y*u.w;
    float   zw = u.z*u.w;
    expr[0] = 
        	  yy*( c1 + c11 + c12 - 2*c14 - 2*c17 + c19 + c2 + c20 - c21 - c23 + c3 + c4 - c6 - c8)
			+ zz*( c0 + c11 - c14 - c16 - c17 + c19 - 2*(c6+c21) + c25 + c26 - c28 + c3 + c4 + c9)
			+ xw*8*( c11 - c17 + c3 - c6)
        +2*(
        	  xx*( 2*c0 - c1 - c10 + c11 + c12 - c13 - c14 + c2 - c3 + c4 + c5 - c6 + c7 - c8 - c9)
			+ ww*( -c0 + 2*c1 - c11 + c13 + c15 - c17 + c19 - c2 - c21 + c25 - c27 - c28 + c3 - c5 + c9)
        )
		+4*(
        	  xy*( c0 + c11 + c12 - c14 - c17 + c2 + c4 - c6 - c8 - c9)
			+ xz*( c0 - c1 + 2*(c11-c6) - c14 - c17 + c2 + c3 + c4 - c9)
			+ yz*( c11 - c14 - c17 + c19 - c21 + c3 + c4 - c6)
			+ yw*( -c0 + c1 + c11 + c19 - c2 - c21 + 2*(c3-c17) - c6 + c9)
			+ zw*( c1 - c17 + c19 - c2 - c21 + c25 - c28 + c3 - c6 + c9)
        );
    expr[1] = 
        	  yy*( c1 - 2*(c12+c19) + c13 + c14 + c17 + c18 + c2 - c21 - c22 + c3 + c4 - c7 - c8)
			+ zz*( -4*(c0+c19-c2) + 4*c1 + 2*c11 + c15 + c16 + c17 + c18 - c21 - c22 - c25 - c26 - c3 - c4 - c5 - c6 + 2*c9)
			+ xw*8*( c13 - c19 + c2 - c7)
        +2*(
        	  ww*( -c0 + 2*c1 + c11 - c13 + c15 + c17 - c19 + c2 - c21 - c25 + c27 - c28 - c3 - c5 + c9)
			+ xx*( 2*c0 - c1 - c10 - c11 - c12 + c13 + c14 - c2 + c3 + c4 + c5 + c6 - c7 - c8 - c9)
			+ yz*( -2*(c0-c9+2*c19) + 3*(c1+c2) + c13 + c14 + c17 + c18 - c21 - c22 - c3 - c4 - c7 - c8 )
        )
		+4*(
        	  xy*( c0 - c12 + c13 + c14 - c19 + c3 + c4 - c7 - c8 - c9)
			+ xz*( c1 + c13 + c14 - 2*c19 + c2 - c7 - c8)
			+ yw*( -c0 + c1 + c13 + c17 - 2*(c19-c2) - c21 - c3 - c7 + c9)
			+ zw*( -2*(c0-c1+c19) + c11 + c15 + c17 + 2*c2 - c21 - c25 - c3 - c5 + c9)
        );
    expr[2] = 
			  yy*( 4*(c0-c1-c4+c9) - c11 - c12 - c13 - c14 - c17 - c18 - c19 + 2*(c2+c3) - c20 + c22 + c23 + c6 + c7 )
			+ zz*( c0 + c11 - c14 - c15 - c18 + c19 + c21 + c22 - 2*(c25+c4) - c28 + c5 + c6 + c9)+
			+ xw*8*( -c15 - c4 + c5 + c9)
        +2*(
        	  xx*( 2*c0 - c1 + c10 - c11 - c12 - c13 - c14 + c2 + c3 - c4 - c5 + c6 + c7 - c8 + c9)
			+ ww*( -c0 + 2*c1 + c11 + c13 - c15 + c17 + c19 - c2 + c21 - c25 - c27 - c28 - c3 + c5 - c9)
			+ yz*( 3*(c0+c9) - 2*(c1-c2+2*c4) - c11 - c14 - c15 - c18 - c19  + c21 + c22 - c28  + c5 + c6 )
        )
		+4*(
        	  xy*( 2*(c0-c1-c4+c9) - c11 - c12 - c13 - c14 + c2 + c3 + c6 + c7 )
			+ xz*( c0 - c1 - c11 - c14 - c15 + c2 - 2*(c4-c9) + c5 + c6 )
			+ yw*( c0 - c15 + c21 - c28 - 2*c4 + c5 + c9)
			+ zw*( c1 + c11 - c15 + c19 - c2 + c21 - c25 - c28 - c4 + c5)
        );
    expr[3] = 
        	  yy*( c1 + c11 + c12 - 2*(c13+c18) + c19 + c2 + c20 - c22 - c24 + c3 + c4 - c5 - c7)
			+ zz*( c0 + c11 - c13 - c15 - c18 + c19 - 2*(c22+c5) + c25 + c26 - c29 + c3 + c4 + c9)
			+ xw*8*( c0 + c1 + c11 - c13 - c2 + c3 - c5 - c9)
        +2*(
        	  xx*( 2*c0 + c1 - c10 + c11 + c12 - c13 - c14 - c2 + c3 - c4 - c5 + c6 - c7 + c8 - c9)
			+ ww*( c0 + 2*c1 + c11 - c13 - c15 + c17 - c19 - c2 - c21 + c25 - c27 + c28 + c3 - c5 - c9)
        )
		+4*(
        	  xy*( c0 + c1 + c11 + c12 - c13 - c18 + c3 - c5 - c7 - c9)
			+ xz*( c0 + c1 + 2*(c11-c5) - c13 - c18 - c2 + c3 + c4 - c9)
			+ yz*( c11 - c13 - c18 + c19 - c22 + c3 + c4 - c5)
			+ yw*( c0 + c1 + c11 - 2*(c13-c3) + c19 - c2 - c22 - c5 - c9)
			+ zw*( c0 + c1 + c11 - c13 - c15 - c2 - c22 + c25 + c3 - c5)
        );
    expr[4] = 
        	  yy*( c1 - 2*(c11+c20) + c13 + c14 + c17 + c18 + c2 - c23 - c24 + c3 + c4 - c5 - c6)
			+ zz*( 2*(c0+c19) + 4*(c1-c11+c2-c9) + c13 + c14  - c21 - c22 - c25 - c26 + c28 + c29 - c3 - c4 - c5 - c6 )
			+ xw*8*( c0 + c1 - c11 + c13 + c2 - c3 - c5 - c9)
		+2*(
        	  xx*( 2*c0 + c1 - c10 - c11 - c12 + c13 + c14 + c2 - c3 - c4 - c5 - c6 + c7 + c8 - c9)
			+ ww*( c0 + 2*c1 - c11 + c13 - c15 - c17 + c19 + c2 - c21 - c25 + c27 + c28 - c3 - c5 - c9)
			+ yz*( 2*(c0-2*c11-c9) + 3*(c1+c2)  + c13 + c14 + c17 + c18  - c23 - c24 - c3 - c4 - c5 - c6 )
        )
		+4*(
        	  xy*( c0 + c1 - c11 + c13 + c14 + c2 - c20 - c5 - c6 - c9)
			+ xz*( 2*(c0+c1-c11+c2-c9) + c13 + c14 - c3 - c4 - c5 - c6 )
			+ yw*( c0 + c1 - 2*(c11-c2) + c13 + c17 - c23 - c3 - c5 - c9)
			+ zw*( c0 + 2*(c1-c11-c9) + c13 + c19 + 2*c2 - c21 - c25 + c28 - c3 - c5 )
        );
    expr[5] = 
        	  yy*( 4*(c0-c2-c3+c9) + 2*(c1+c4) - c11 - c12 - c13 - c14 - c17 - c18 - c19  - c20 + c21 + c24  + c5 + c8)
			+ zz*( c0 + c11 - c13 - c16 - c17 + c19 + c21 + c22 - 2*c26 - c29 - 2*c3 + c5 + c6 + c9)+
			+ xw*8*( c0 + c1 - c11 - c13 - c2 - c3 + c5 + c9)
        +2*(
        	  xx*( 2*c0 + c1 + c10 - c11 - c12 - c13 - c14 - c2 - c3 + c4 + c5 - c6 - c7 + c8 + c9)
			+ ww*( c0 + 2*c1 - c11 - c13 + c15 - c17 - c19 - c2 + c21 - c25 - c27 + c28 - c3 + c5 + c9)
			+ yz*( 3*(c0+c9) + 2*(c1-c2) - c11 - c13 - c16 - c17 - c19  + c21 + c22 - c29 - 4*c3 + c5 + c6 )
        )
		+4*(
        	  xy*( 2*(c0-c2-c3+c9) + c1 - c11 - c12 - c13 - c14 + c4 + c5 + c8 )
			+ xz*( c0 + c1 - c11 - c13 - c16 - c2 - 2*(c3-c9) + c5 + c6 )
			+ yw*( 2*(c0+c1-c2-c3+c9) - c11 - c13 - c17 - c19 + c21 + c5 )
			+ zw*( c0 + c1 - c13 - c17 - c2 + c21 - c26 - c3 + c5 + c9)
        );
    for(int i=0 ; i<6 ; i++)    expr[i] *= SCALE_Tij;
    return expr;
}

const int map_Tij[24*6] = 
int[]
(
    4, 3, 2, 5, 0, 1,   //  0
    3, 4, 2, 5, 1, 0,   //  1
    4, 5, 0, 3, 2, 1,   //  2
    3, 5, 1, 4, 2, 0,   //  3
    5, 3, 1, 4, 0, 2,   //  4
    5, 4, 0, 3, 1, 2,   //  5

    1, 0, 2, 5, 3, 4,   //  6
    0, 1, 2, 5, 4, 3,   //  7
    1, 2, 0, 3, 5, 4,   //  8
    0, 2, 1, 4, 5, 3,   //  9
    2, 0, 1, 4, 3, 5,   //  10
    2, 1, 0, 3, 4, 5,   //  11

    1, 3, 5, 2, 0, 4,   //  12
    0, 4, 5, 2, 1, 3,   //  13
    1, 5, 3, 0, 2, 4,   //  14
    0, 5, 4, 1, 2, 3,   //  15
    2, 3, 4, 1, 0, 5,   //  16
    2, 4, 3, 0, 1, 5,   //  17

    4, 0, 5, 2, 3, 1,   //  18
    3, 1, 5, 2, 4, 0,   //  19
    4, 2, 3, 0, 5, 1,   //  20
    3, 2, 4, 1, 5, 0,   //  21
    5, 0, 4, 1, 3, 2,   //  22
    5, 1, 3, 0, 4, 2    //  23
);

float[6] compute_Hessian(vec3 p)
{
    float    expr[6] = eval_Tij();

    float	d12 = expr[map_Tij[6*idx + 0]];
    float	d13 = expr[map_Tij[6*idx + 1]];
    float	d14 = expr[map_Tij[6*idx + 2]];
    float	d23 = expr[map_Tij[6*idx + 3]];
    float	d24 = expr[map_Tij[6*idx + 4]];
    float	d34 = expr[map_Tij[6*idx + 5]];

    float	Dxx = 0.25*(-d12-d13        -d24-d34)*scale_norm.x*scale_norm.x;
    float	Dyy = 0.25*(-d12    -d14-d23    -d34)*scale_norm.y*scale_norm.y;
    float	Dzz = 0.25*(    -d13-d14-d23-d24    )*scale_norm.z*scale_norm.z;
    float	Dyz = 0.25*(        -d14+d23        )*scale_norm.y*scale_norm.z;
    float	Dzx = 0.25*(     d13        -d24    )*scale_norm.x*scale_norm.z;
    float	Dxy = 0.25*( d12                -d34)*scale_norm.x*scale_norm.y;
    return float[6](Dxx,Dyy,Dzz,Dyz,Dzx,Dxy);
}

#define	SHADING_BLINN_PHONG 1
struct TMaterial
{
	vec3	ambient;
	vec3	diffuse;
	vec3	specular;
	vec3	emission;
	float	shininess;
};
struct TLight
{
	vec4	position;
	vec3	ambient;
	vec3	diffuse;
	vec3	specular;
};

TLight		uLight = TLight(
        vec4(1,1,1,0),
        vec3(.2,.2,.2),
        vec3(1,1,1),
        vec3(1,1,1)
        );

vec4 shade_Blinn_Phong(vec3 n, vec4 pos_eye, TMaterial material, TLight light)
{
	vec3	l;
	if(light.position.w == 1.0)
		l = normalize((light.position - pos_eye).xyz);		// positional light
	else
		l = normalize((light.position).xyz);	// directional light
	vec3	v = -normalize(pos_eye.xyz);
	vec3	h = normalize(l + v);
	float	l_dot_n = max(dot(l, n), 0.0);
	vec3	ambient = light.ambient * material.ambient;
	vec3	diffuse = light.diffuse * material.diffuse * l_dot_n;
	vec3	specular = vec3(0.0);

	if(l_dot_n >= 0.0)
	{
		specular = light.specular * material.specular * pow(max(dot(h, n), 0.0), material.shininess);
	}
	return vec4(ambient + diffuse + specular, 1);
}


vec4 compute_color(vec4 p, vec3 g, float[6] d2) {

	float	Dxx = d2[0];
	float	Dyy = d2[1];
	float	Dzz = d2[2];
	float	Dyz = d2[3];
	float	Dzx = d2[4];
	float	Dxy = d2[5];

	mat3	H = mat3(Dxx, Dxy, Dzx,
					Dxy, Dyy, Dyz,
					Dzx, Dyz, Dzz);
	float	one_over_len_g = 1.0/length(g);
	vec3	n = -g*one_over_len_g;
	mat3    P = mat3(1.0) - mat3(n.x*n.x, n.x*n.y, n.x*n.z,
								n.x*n.y, n.y*n.y, n.y*n.z,
								n.x*n.z, n.y*n.z, n.z*n.z);
	mat3    M = -P*H*P*one_over_len_g;
	float   T = M[0][0] + M[1][1] + M[2][2];
	mat3    MMt = M*transpose(M);
	float   F = sqrt(MMt[0][0] + MMt[1][1] + MMt[2][2]);
	float   k_max = (T + sqrt(2.0*F*F - T*T))*0.5;
	float   k_min = (T - sqrt(2.0*F*F - T*T))*0.5;

	float	scale_k = 10;
	vec2	tc = vec2(scale_k*vec2(k_max,k_min)+0.5);

	if(p.w!=0.0)
	{
#if	SHADING_BLINN_PHONG
		TMaterial	material = 
			TMaterial(
				vec3(.1,.1,.1),
				texture(tex_colormap_2d, tc).xyz,
				vec3(1,1,1),
				vec3(0,0,0),
				128.0*0.5
				);
		return shade_Blinn_Phong(normalize(mat3(MV)*(-p.w*g)), MV*vec4(p.xyz,1), material, uLight);
#else
		return texture(tex_colormap_2d, tc);
#endif
	}
}



void main() {

	vec3 start = texture(tex_front, vTexCoord).xyz*scale_lattice + offset_lattice;
	vec3 end = texture(tex_back, vTexCoord).xyz*scale_lattice + offset_lattice;

	vec3	p = start;
	vec3	p_prev;
	vec3	dir = normalize(end-start);

	float	step = scale_step*dim_max;

	float	len = 0;
	float	len_full = length(end - start);
	float	voxel, voxel_prev;

    voxel = EVAL(p);

    float   orientation = 2.0*float(voxel < level)-1.0;	// equivalent to (voxel<level?1:-1)

	for(int i = 0 ; i < 1000 ; i++)
	{
		p += step*dir;
		len += step;
		if(len > len_full)
		{
			DISCARD;
		}

		voxel = EVAL(p);

		if(orientation*voxel > orientation*level)
		{
			// One step of Regula Falsi
			if(abs(voxel-voxel_prev) > 0.00001)
			{
				p = (p*(voxel_prev-level) - p_prev*(voxel-level))/(voxel_prev-voxel);
#ifndef NO_PREPROCESS
				preprocess(p);
#endif
#ifndef NO_FETCH_COEFFICIENTS
				fetch_coefficients();
#endif
			}
            vec3    g = compute_gradient(p);
            float[6] H = compute_Hessian(p);

            vec3        pos = (p*vec3(scale_norm) + vec3(offset_norm) - vec3(.5,.5,.5));
            fColor = compute_color(vec4(pos,orientation), g, H);
            return;
 
	}
		voxel_prev = voxel;
		p_prev = p;
	}
	DISCARD;
}

