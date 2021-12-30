#version 3.5;

#declare param=array[3];
#fopen File "./dim.dat" read
#declare i=0;
#while (defined(File))
    #read (File,l)
    #declare param[i]=l;
    #declare i=i+1;
#end
#fclose File

#declare nx=param[0];
#declare ny=param[1];
#declare nz=param[2];
#declare Np=nx*ny*nz;
#declare spins=array[Np];

#declare cama=array[2];
#fopen File "./cam_temp.dat" read
#declare i=0;
#while (defined(File))
    #read (File,l)
    #declare cama[i]=l;
    #declare i=i+1;
#end
#fclose File


#declare r=1.5*sqrt(nx*nx+ny*ny+nz*nz);
#declare phi=pi/6;
#declare cx=r*cama[0]*cos(phi);
#declare cy=r*cama[1]*cos(phi);
#declare cz=r*sin(phi);
#declare cam=<cx,cy,cz>;
#declare midspace=<nx/2,ny/2,nz/2>;
#declare cam=midspace+cam;

#declare i=0;
#while(i<Np)
    #declare spins[i]=0;
    #declare i=i+1;
#end
//-------------
background
{
    color rgb<0,0,0>
}
camera
{
    location cam
    sky<0,0,1>
    right<-1,0,0>
    look_at midspace
}
light_source
{
    cam
    color rgb<1,1,1>
    spotlight radius 1
    adaptive 1
    jitter
    point_at midspace
    shadowless
}
//--------------


#declare tr=0.7;
#declare up_top=
difference
{
    sphere{<0,0,0>,0.5}
    box{<-0.5,-0.5,-0.5>,<0.5,0.5,0>}
    texture { pigment { color rgb <1,0,0> transmit tr} }
    finish  { phong 1.0 phong_size 60}
}
#declare up_bottom=
difference
{
    sphere{<0,0,0>,0.5}
    box{<-0.5,-0.5,0>,<0.5,0.5,0.5>}
    texture { pigment { color rgb <1,1,1> transmit tr} }
    finish  { phong 1.0 phong_size 60}
}
#declare up_spin=
union
{
    object{up_top}
    object{up_bottom}
}
//----------------
#declare down_top=
difference
{
    sphere{<0,0,0>,0.5}
    box{<-0.5,-0.5,-0.5>,<0.5,0.5,0>}
    texture { pigment { color rgb <1,1,1> transmit tr}}
    finish  { phong 1.0 phong_size 60 }
}
#declare down_bottom=
difference
{
    sphere{<0,0,0>,0.5}
    box{<-0.5,-0.5,0>,<0.5,0.5,0.5>}
    texture { pigment { color rgb <1,0,0> transmit tr}}
    finish  { phong 1.0 phong_size 60}
}
#declare down_spin=
union
{
    object{down_top}
    object{down_bottom}
}
//----------------

#fopen MyFile "./frame_data.dat" read
#declare i=0;
#while (defined(MyFile))
    #read (MyFile,var1)
    #declare spins[i]=var1;
    #declare i=i+1;
#end
#fclose MyFile

#declare k=0;
#while(k<nz)
    #declare j=0;
    #declare kk=k*nx*ny;
    #while(j<ny)
        #declare i=0;
        #declare jj=j*nx;
        #while(i<nx)
            #declare ii=kk+jj+i;
            #if(spins[ii]=1)
                object{up_spin translate<i,j,k>}
            #else
                #if(spins[ii]=-1)
                    object{down_spin translate<i,j,k>}
                #end
            #end
            #declare i=i+1;
        #end
        #declare j=j+1;
    #end
    #declare k=k+1;
#end


