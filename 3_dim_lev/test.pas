program ex2dimarray;

uses
  sysutils;


var 
   a: array [0..3, 0..3] of integer;
   i,j,k : integer;

begin
   for i:=0 to 3 do
      for j:=0 to 3 do
         a[i,j]:= i * j;  
   
   for k:=0 to 1000 do   
      for i:=0 to 3 do
      begin  
         for j:=0 to 3 do  
            write(a[i,j]:2,' ');  
         writeln;  
      end;
   write(maxint);
   writeln;  
end.