{ (c) W. J. Heeringa 2021 }

program calculate_word_distance(input,output,stderr);

uses
  sysutils;
  
const
  ml = 1000;
  mp = 1000;
  mt = 1000;

type
  stresstype     = array[1..2] of char;
  headtype       = array[1..3] of char;
  supratype      = array[1..2] of char;
  diacritictype  = array[1..2] of char;

  segmenttype    = record
                     head      : headtype;
                     supra     : array[1..4] of supratype;
                     diacritic : array[1..4] of diacritictype;
                   end;

  lettertype     = record
                     stress    : stresstype;
                     segment   : array[1..2] of segmenttype;
                   end;

  wordtype       = record
                      letter    : array[1..80] of lettertype;
                     number    : integer;
                   end;

  comptype      = record
                    number    : integer;
                    segment   : array[1..mp] of segmenttype;
                    matrix    : array[1..mp] of array[1..mp] of real;
                  end;

  segmenttype0  = record
                    sound     : integer;
                  end;
                
  lettertype0   = record
                    segment   : array[1..2] of segmenttype0;
                  end;
                
  wordtype0     = record
                    letter    : array[1..80] of lettertype0;
                    number    : integer;
                  end;

  matrixtype     = array[0..80,0..80] of real;

  tracetype      = array[0..80,0..80] of integer;

  alignrec       = record
                     word1     :  integer;
                     word2     :  integer;
                     dist      :  real;
                   end;

  aligntype      = array[1..80] of alignrec;

var
  comp               : comptype;
  vc                 : real;
  empty              : lettertype0;
  matrix1,matrix2    : matrixtype;
  trace              : tracetype;
  align              : aligntype;
  stringO,stringY    : string[255];
  wordO  ,wordY      : wordtype;
  wordO0 ,wordY0     : wordtype0;
  cfile              : string[255];
  part               : integer;
  dist0, dist1, dist2: real;

// Functions for showing and getting parameters

procedure usage;
begin
  writeln(stderr);
  writeln(stderr,'(c) W. J. Heeringa 2021');
  writeln(stderr);
  writeln(stderr,'Usage: ',ParamStr(0),' dialectA dialectB cmpfile 1-7');
  writeln(stderr);
  writeln(stderr,'dialectA: X-SAMPA transcription of a word in dialect A');
  writeln(stderr,'dialectB: X-SAMPA transcription of a word in dialect B');
  writeln(stderr);
  writeln(stderr,'cmpfile : file with segment comparisons'); 
  writeln(stderr);
  writeln(stderr,'1: Only vowels                               ');
  writeln(stderr,'2: Only vowel     vs. vowel     substitutions');
  writeln(stderr,'3: Only vowel     indels                     ');
  writeln(stderr,'4: Only consonants                           ');
  writeln(stderr,'5: Only consonant vs. consonant substitutions');
  writeln(stderr,'6: Only consonant indels                     ');
  writeln(stderr,'7: Only schwa     vs. sonorant  substitutions');
  writeln(stderr);
 
  writeln(stderr);
  halt;
end;

procedure getparameters;
begin
  if (paramcount=3) or (paramcount=4)
    then begin
           stringO:=paramstr(1);
           stringY:=paramstr(2);
           cfile  :=paramstr(3);
         end
    else usage;
    
  if (paramcount=3)
    then part:=0
    else begin     
           if paramstr(4)='1' then part:=1 else
           if paramstr(4)='2' then part:=2 else
           if paramstr(4)='3' then part:=3 else
           if paramstr(4)='4' then part:=4 else
           if paramstr(4)='5' then part:=5 else
           if paramstr(4)='6' then part:=6 else
           if paramstr(4)='7' then part:=7 else usage;
         end;
end;

// Functions for parsing input transcriptions

function instresss(ch:char):boolean;
begin
  if ch=''''
    then instresss:=true
    else
  if ch='"'
    then instresss:=true
    else
  if ch='%'
    then instresss:=true
    else instresss:=false;
end;

function inheads(ch:char):boolean;
begin
  if ch='0'
    then inheads:=true
    else
  if ch='i'
    then inheads:=true
    else
  if ch='y'
    then inheads:=true
    else
  if ch='1'
    then inheads:=true
    else
  if ch='}'
    then inheads:=true
    else
  if ch='M'
    then inheads:=true
    else
  if ch='u'
    then inheads:=true
    else
  if ch='I'
    then inheads:=true
    else
  if ch='Y'
    then inheads:=true
    else
  if ch='U'
    then inheads:=true
    else
  if ch='e'
    then inheads:=true
    else
  if ch='2'
    then inheads:=true
    else
  if ch='8'
    then inheads:=true
    else
  if ch='7'
    then inheads:=true
    else
  if ch='o'
    then inheads:=true
    else
  if ch='@'
    then inheads:=true
    else
  if ch='E'
    then inheads:=true
    else
  if ch='9'
    then inheads:=true
    else
  if ch='3'
    then inheads:=true
    else
  if ch='V'
    then inheads:=true
    else
  if ch='O'
    then inheads:=true
    else
  if ch='{'
    then inheads:=true
    else
  if ch='6'
    then inheads:=true
    else
  if ch='a'
    then inheads:=true
    else
  if ch='&'
    then inheads:=true
    else
  if ch='A'
    then inheads:=true
    else
  if ch='Q'
    then inheads:=true
    else
  if ch='p'
    then inheads:=true
    else
  if ch='b'
    then inheads:=true
    else
  if ch='t'
    then inheads:=true
    else
  if ch='d'
    then inheads:=true
    else
  if ch='c'
    then inheads:=true
    else
  if ch='k'
    then inheads:=true
    else
  if ch='g'
    then inheads:=true
    else
  if ch='q'
    then inheads:=true
    else
  if ch='?'
    then inheads:=true
    else
  if ch='m'
    then inheads:=true
    else
  if ch='F'
    then inheads:=true
    else
  if ch='n'
    then inheads:=true
    else
  if ch='J'
    then inheads:=true
    else
  if ch='N'
    then inheads:=true
    else
  if ch='r'
    then inheads:=true
    else
  if ch='4'
    then inheads:=true
    else
  if ch='B'
    then inheads:=true
    else
  if ch='f'
    then inheads:=true
    else
  if ch='v'
    then inheads:=true
    else
  if ch='T'
    then inheads:=true
    else
  if ch='D'
    then inheads:=true
    else
  if ch='s'
    then inheads:=true
    else
  if ch='z'
    then inheads:=true
    else
  if ch='S'
    then inheads:=true
    else
  if ch='Z'
    then inheads:=true
    else
  if ch='C'
    then inheads:=true
    else
  if ch='x'
    then inheads:=true
    else
  if ch='G'
    then inheads:=true
    else
  if ch='X'
    then inheads:=true
    else
  if ch='R'
    then inheads:=true
    else
  if ch='h'
    then inheads:=true
    else
  if ch='K'
    then inheads:=true
    else
  if ch='w'
    then inheads:=true
    else
  if ch='P'
    then inheads:=true
    else
  if ch='j'
    then inheads:=true
    else
  if ch='M'
    then inheads:=true
    else
  if ch='l'
    then inheads:=true
    else
  if ch='L'
    then inheads:=true
    else inheads:=false;
end;

function insupras(ch:char):boolean;
begin
  if ch=':'
    then insupras:=true
    else insupras:=false;
end;

function indiacritics(ch:char):boolean;
begin
  if ch='0'
    then indiacritics:=true
    else
  if ch='v'
    then indiacritics:=true
    else
  if ch='h'
    then indiacritics:=true
    else
  if ch='+'
    then indiacritics:=true
    else
  if ch='-'
    then indiacritics:=true
    else
  if ch='='
    then indiacritics:=true
    else
  if ch='t'
    then indiacritics:=true
    else
  if ch='k'
    then indiacritics:=true
    else
  if ch='w'
    then indiacritics:=true
    else
  if ch='j'
    then indiacritics:=true
    else
  if ch='G'
    then indiacritics:=true
    else
  if ch='?'
    then indiacritics:=true
    else
  if ch='e'
    then indiacritics:=true
    else
  if ch='a'
    then indiacritics:=true
    else
  if ch='~'
    then indiacritics:=true
    else
  if ch='}'
    then indiacritics:=true
    else
  if ch='r'
    then indiacritics:=true
    else
  if ch='o'
    then indiacritics:=true
    else
  if ch='_'
    then indiacritics:=true
    else
  if ch='X'
    then indiacritics:=true
    else indiacritics:=false;
end;

procedure halt0(mess:string);
begin
  writeln(stderr,mess);
  halt;
end;

procedure initstress(var stress:stresstype);
begin
  stress:='  ';
end;

procedure initsegment(var segment:segmenttype);
var
  p,d : integer;
begin
  segment.head:='   ';

  for p:=1 to 4 do begin
    segment.supra[p]:='  ';
  end;

  for d:=1 to 4 do begin
    segment.diacritic[d]:='  ';
  end;
end;

procedure initletter(var letter:lettertype);
var
  s : integer;
begin
  initstress(letter.stress);

  for s:=1 to 2 do begin
    initsegment(letter.segment[s]);
  end;
end;

procedure rstress(stringX:string;var stress:stresstype;var pos:integer);
begin
  if (pos<=Length(stringX))
    then begin
           if stringX[pos]=''''
             then begin
                    stress[1]:=stringX[pos];
                    inc(pos);
                  end
             else
           if stringX[pos]='"'
             then begin
                    stress[1]:=stringX[pos];
                    inc(pos);
                    if (pos<=Length(stringX)) and (stringX[pos]='"')
                      then begin
                             stress[2]:=stringX[pos];
                             inc(pos)
                           end
                      else {nothing}
                  end
             else
           if stringX[pos]='%'
             then begin
                    stress[1]:=stringX[pos];
                    inc(pos);
                    if (pos<=Length(stringX)) and (stringX[pos]='%')
                      then begin
                             stress[2]:=stringX[pos];
                             inc(pos);
                           end
                      else {nothing}
                  end
             else halt0(Concat(stringX,' ',IntToStr(pos),' '' or % or " or ^ expected'));
         end
    else {nothing};              
end;

procedure rhead(stringX:string;var head:headtype;var pos:integer);
begin
  if inheads(stringX[pos])
    then begin
           head[1]:=stringX[pos];
            
           inc(pos);
           if (pos<=Length(stringX)) and (stringX[pos]='\')
             then begin
                    head[2]:=stringX[pos];
                    inc(pos);
                  end
             else {nothing};
           if (pos<=Length(stringX)) and (stringX[pos]='`')
             then begin
                    head[3]:=stringX[pos];
                    inc(pos);
                  end
             else {nothing}
         end
    else halt0(Concat(stringX,' ',IntToStr(pos),' head expected'))
end;

procedure rsupra(stringX:string;var supra:supratype;var pos:integer);
begin
  if stringX[pos]=':'
    then begin
           supra[1]:=stringX[pos];
           inc(pos);
         end
    else halt0(Concat(stringX,' ',IntToStr(pos),' : expected'));

  if stringX[pos]='\'
    then begin
           supra[2]:=stringX[pos];
           inc(pos);
         end
    else {nothing};
end;

procedure rdiacritic(stringX:string;var diacritic:diacritictype;var pos:integer);
begin
  if stringX[pos]='_'
    then begin
           inc(pos);
           
           if (pos<=Length(stringX)) and (indiacritics(stringX[pos]))
             then begin
                    diacritic[1]:=stringX[pos];
                    inc(pos);
                  end
             else halt0(Concat(stringX,' ',IntToStr(pos),' diacritic expected'));

           if (pos<=Length(stringX)) and (stringX[pos]='\')
             then begin
                    diacritic[2]:=stringX[pos];
                    inc(pos);
                  end
             else {nothing};
         end
    else halt0(Concat(stringX,' ',IntToStr(pos),' _ expected'));
end;

procedure rsegment(stringX:string;var segment:segmenttype;var pos:integer);
var
  p,d : integer;
begin
  p:=0;
  d:=0;

  rhead(stringX,segment.head,pos);
  
  if (pos<=Length(stringX))
    then while insupras(stringX[pos]) or (stringX[pos]='_') do begin
           if insupras(stringX[pos])
             then begin
                    inc(p);
                    rsupra(stringX,segment.supra[p],pos)
                  end
             else begin
                    inc(d);
                    rdiacritic(stringX,segment.diacritic[d],pos);
                  end;
         end
    else {nothing}; 
end;

procedure rletter(stringX:string;var letter:lettertype;var pos:integer);
var
  s : integer;
begin
  s:=1;

  initletter(letter);
  if instresss(stringX[pos]) or inheads(stringX[pos])
    then begin
           if instresss(stringX[pos])
             then rstress(stringX,letter.stress,pos)
             else {nothing};
            
           rsegment(stringX,letter.segment[s],pos)
         end
    else if stringX[pos]='['
          then begin
                 inc(pos);
                 
                 if (pos<=Length(stringX))
                   then begin
                          if instresss(stringX[pos])
                            then rstress(stringX,letter.stress,pos)
                            else {nothing};
                          
                          if (pos<=Length(stringX))
                            then rsegment(stringX,letter.segment[s],pos)
                            else {nothing};
                          
                          if (pos<=Length(stringX))
                            then begin
                                   while (pos<=Length(stringX)) and inheads(stringX[pos]) do begin
                                     inc(s);
                                     rsegment(stringX,letter.segment[s],pos)
                                   end;

                                   if stringX[pos]=']'
                                     then inc(pos)
                                     else halt0(Concat(stringX,' ',IntToStr(pos),' ] expected'));
                                 end
                            else {nothing}; 
                        end
                   else {nothing};                        
               end
          else halt0(Concat(stringX,' ',IntToStr(pos),' stress or head or [ expected'));
end;

procedure rword(stringX:string;var w_rd:wordtype);
var
  l, pos : integer;
begin
  l  := 1;
  pos:= 1;
  
  if (pos<=Length(stringX))
    then begin
           rletter(stringX,w_rd.letter[l],pos);
           while (pos<=Length(stringX)) and (instresss(stringX[pos]) or inheads(stringX[pos]) or (stringX[pos]='[')) do begin
             inc(l);
             rletter(stringX,w_rd.letter[l],pos);
           end;
           w_rd.number:=l;
         end
    else w_rd.number:=0;   
end;

// Functions for reading segment distances and finding the maximum

procedure openr(var fp:text;name:string);
begin
  assign(fp,name);
  {$I-}
  reset(fp);
  {$I+}

  if IOResult <> 0
    then begin
           writeln(stderr,'Error opening file ',name);
           halt;
         end
    else {nothing};
end;

procedure initcomp;
var
  fp      : text;
  i,j,pos : integer;
  v       : real;
  stringX : string[255];
begin
  openr (fp,cfile);
  readln(fp,comp.number);

  for i:=1 to comp.number do begin
    initsegment(comp.segment[i]);
    readln(fp,stringX);
    pos:=1;
    rsegment(stringX,comp.segment[i],pos);
    comp.matrix[i,i]:=0;
  end;

  for i:=2 to comp.number do begin
    for j:=1 to (i-1) do begin
      readln(fp,v);
      comp.matrix[i,j]:=v;
      comp.matrix[j,i]:=v;
    end;
  end;

  close(fp);
end;

procedure initvc;
var
  i,j : integer;
begin
  vc:=0;
  for i:=2 to comp.number do begin
    for j:=1 to (i-1) do begin
      if comp.matrix[i,j]>vc
        then vc:=comp.matrix[i,j]
        else {nothing};
    end;
  end;
end;

// Functions for calculating the distances

function sound(segment:segmenttype):integer;
var
  t     : integer;
  found : boolean;
begin
  t:=0;
  found:=false;

  while (not found) and (t<comp.number) do begin
    t:=t+1;
    found:=(segment.head=comp.segment[t].head);
  end;

  if not found
    then begin
           writeln(stderr,segment.head,' not found in comparison table');
           halt;
         end
    else {nothing};

  sound:=t;
end;

procedure change(w_rd:wordtype;var word0:wordtype0);
var
  l : integer;
begin
  for l:=1 to w_rd.number do begin
    word0.letter[l].segment[1].sound:=sound(w_rd.letter[l].segment[1]);

    if w_rd.letter[l].segment[2].head='   '
      then word0.letter[l].segment[2].sound:=word0.letter[l].segment[1].sound
      else word0.letter[l].segment[2].sound:=sound(w_rd.letter[l].segment[2]);
  end;

  word0.number:=w_rd.number;
end;

procedure initempty;
var
  s : integer;
begin
  for s:=1 to 2 do begin
    empty.segment[s].sound:=1;
  end;
end;

function vowel1(sound:integer):boolean;
begin
  vowel1:=
     ((sound>=02) and (sound<=29));
end;

function vowel10(sound:integer):boolean;
begin
  vowel10:=
     ((sound>=02) and (sound<=29))
end;

function consonant1(sound:integer):boolean;
begin
  consonant1:=
     ((sound>=30) and (sound<=88));
end;

function consonant10(sound:integer):boolean;
begin
  consonant10:=
     ((sound>=30) and (sound<=88))
end;

function vowel2(sound:integer):boolean;
begin
  vowel2:=
      (sound =17);
end;

function consonant2(sound:integer):boolean;
begin
  consonant2:=
     ((sound>=43) and (sound<=49)) or
     ((sound>=50) and (sound<=54)) or
     ((sound>=79) and (sound<=88))
end;

function empty0(sound:integer):boolean;
begin
  empty0:=
      (sound = 1);
end;

function weight0(segment1,segment2:segmenttype0):real;
begin
  if (    vowel1(segment1.sound) and     vowel1(segment2.sound)) or
     (consonant1(segment1.sound) and consonant1(segment2.sound)) or
     (    vowel1(segment1.sound) and     empty0(segment2.sound)) or
     (    empty0(segment1.sound) and     vowel1(segment2.sound)) or
     (consonant1(segment1.sound) and     empty0(segment2.sound)) or
     (    empty0(segment1.sound) and consonant1(segment2.sound)) or
     (    empty0(segment1.sound) and     empty0(segment2.sound))

    then weight0:=comp.matrix[segment1.sound,segment2.sound]
    else weight0:=maxint;
end;

function weight(letter1,letter2:lettertype0):real;
begin
  weight:=(weight0(letter1.segment[1],letter2.segment[1])+
           weight0(letter1.segment[2],letter2.segment[2]))/2;
end;

function min(a,b,c:real):real;
begin
  if ((a<=b) and (a<=c)) then min:=a else
  if ((b<=a) and (b<=c)) then min:=b else min:=c;
end;

function max(a,b,c:real):real;
begin
  if ((a>=b) and (a>=c)) then max:=a else
  if ((b>=a) and (b>=c)) then max:=b else max:=c;
end;

procedure check(var value1:real;prevA1,prevB1,prevC1:real;
                var value2:real;prevA2,prevB2,prevC2:real);
begin
  value1:=min(prevA1,prevB1,prevC1);

  if prevA1<>value1
    then prevA2:=-maxint;

  if prevB1<>value1
    then prevB2:=-maxint;

  if prevC1<>value1
    then prevC2:=-maxint;

  value2:=max(prevA2,prevB2,prevC2);
end;

procedure addtrace(prevA1,prevB1,prevC1,
                   prevA2,prevB2,prevC2:real;l1,l2:integer);
var
  pointer : integer;
begin
  pointer:=0;

  // if (prevA1=matrix1[l1,l2]) and (prevA2=matrix2[l1,l2])
  //   then pointer:=pointer+8;

  // if (prevB1=matrix1[l1,l2]) and (prevB2=matrix2[l1,l2])
  //   then pointer:=pointer+4;

  // if (prevC1=matrix1[l1,l2]) and (prevC2=matrix2[l1,l2])
  //   then pointer:=pointer+2;

  if (prevA1=matrix1[l1,l2])
    then pointer:=pointer+8;

  if (prevB1=matrix1[l1,l2])
    then pointer:=pointer+4;

  if (prevC1=matrix1[l1,l2])
    then pointer:=pointer+2;

  trace[l1,l2]:=pointer;
end;

procedure addalign(l1,l2:integer;d:real;a:integer);
begin
  align[a].word1:=l1;
  align[a].word2:=l2;
  align[a].dist :=d ;
end;

function segment(sound:integer):string;
var
  t     : integer;
  found : boolean;
begin
  t:=0;
  found:=false;

  while (not found) and (t<comp.number) do begin
    t:=t+1;
    found:=(sound=t);
  end;

  if not found
    then begin
           writeln(stderr,sound,' not found in comparison table');
           halt;
         end
    else {nothing};
    
  segment:=comp.segment[t].head;
end;

procedure printalign(a:integer; word1,word2:wordtype0);
var
  l : integer;
begin
  write('*');
  // writeln(stdout);
  for l:=(a-1) downto 1 do begin
    if align[l].word1=0 
      then write('-  ':6)
      else write(segment(word1.letter[align[l].word1].segment[1].sound):6)
  end;

  // writeln(stdout);
  write(#13#10);

  for l:=(a-1) downto 1 do begin
    if align[l].word2=0 
      then write('-  ':6)
      else write(segment(word2.letter[align[l].word2].segment[1].sound):6)
  end;

  write(#13#10);

  // writeln(stdout);
  // // write(stdout, ' ');

  // for l:=(a-1) downto 1 do begin
  //   write(align[l].dist:6:2);
  // end;
  
  // writeln(stdout);
  // writeln(stdout);
end;

procedure procalign(a:integer;word1,word2:wordtype0;var distance1,distance2:real);
var
  l,l0,l1,l2 : integer;
  distance10 : real;
begin
  l :=a;
  l0:=0;
  distance10:=0;

  while (l>1) do begin

    l :=l -1;
    l0:=l0+1;

    if align[l].word1=0 
      then l1:=1
      else l1:=word1.letter[align[l].word1].segment[1].sound;

    if align[l].word2=0 
      then l2:=1 
      else l2:=word2.letter[align[l].word2].segment[1].sound;

    if part=0
      then distance10:=distance10+align[l].dist
      else

    if part=1
      then if (vowel1(l1) and vowel1(l2)) and not(consonant10(l1) and consonant10(l2))
             then distance10:=distance10+align[l].dist
             else
           if ((l1=1) and vowel1(l2)) and not(consonant10(l2))
             then distance10:=distance10+align[l].dist
             else
           if (vowel1(l1) and (l2=1)) and not(consonant10(l1))
             then distance10:=distance10+align[l].dist
             else {nothing}
      else

    if part=2
      then if (vowel1(l1) and vowel1(l2)) and not(consonant10(l1) and consonant10(l2))
             then distance10:=distance10+align[l].dist
             else {nothing}
      else

    if part=3
      then if ((l1=1) and vowel1(l2)) and not(consonant10(l2))
             then distance10:=distance10+align[l].dist
             else
           if (vowel1(l1) and (l2=1)) and not(consonant10(l1))
             then distance10:=distance10+align[l].dist
             else {nothing}
      else

    if part=4
      then if (consonant1(l1) and consonant1(l2)) and not(vowel10(l1) and vowel10(l2))
             then distance10:=distance10+align[l].dist
             else
           if ((l1=1) and consonant1(l2)) and not(vowel10(l2))
             then distance10:=distance10+align[l].dist
             else
           if (consonant1(l1) and (l2=1)) and not(vowel10(l1))  
             then distance10:=distance10+align[l].dist
             else {nothing}
      else

    if part=5
      then if (consonant1(l1) and consonant1(l2)) and not(vowel10(l1) and vowel10(l2))
             then distance10:=distance10+align[l].dist
             else {nothing}
      else

    if part=6
      then if ((l1=1) and consonant1(l2)) and not(vowel10(l2))
             then distance10:=distance10+align[l].dist
             else
           if (consonant1(l1) and (l2=1)) and not(vowel10(l1))  
             then distance10:=distance10+align[l].dist
             else {nothing}
      else

    if part=7
      then if (vowel2(l1) and consonant2(l2)) and not(vowel10(l2))
             then distance10:=distance10+align[l].dist
             else
           if (consonant2(l1) and vowel2(l2)) and not(vowel10(l1))
             then distance10:=distance10+align[l].dist
             else {nothing}
      else usage;
  end;

  distance1:=distance10/vc;
  distance2:=l0;
  
  printalign(a, word1, word2);
  // write('- -');
end;

procedure proctrace(l1,l2,a:integer;word1,word2:wordtype0;var distance1,distance2:real);
var
  pointer : integer;
begin
  if (l1>0) or (l2>0)
    then begin
           pointer:=trace[l1,l2];

           if ((pointer= 8) or (pointer=10) or (pointer=12) or (pointer=14)) and (l1>0)
             then begin
                    addalign(l1,0,matrix1[l1,l2]-matrix1[l1-1,l2],a);
                    proctrace(l1-1,l2,a+1,word1,word2,distance1,distance2);
                    pointer:=pointer-8;
                  end;

           if ((pointer= 4) or (pointer= 6) or (pointer=12) or (pointer=14)) and (l1>0) and (l2>0)
             then begin
                    addalign(l1,l2,matrix1[l1,l2]-matrix1[l1-1,l2-1],a);
                    proctrace(l1-1,l2-1,a+1,word1,word2,distance1,distance2);
                    pointer:=pointer-4;
                  end;

           if ((pointer= 2) or (pointer= 6) or (pointer=10) or (pointer=14)) and (l2>0)
             then begin
                    addalign(0,l2,matrix1[l1,l2]-matrix1[l1,l2-1],a);
                    proctrace(l1,l2-1,a+1,word1,word2,distance1,distance2);
                    pointer:=pointer-2;
                  end;
         end
    else procalign(a,word1,word2,distance1,distance2);
end;

procedure determine(word1,word2:wordtype0;var distance0, distance1,distance2:real);
var
  l1,l2                : integer;
  prevA1,prevB1,prevC1 : real;
  prevA2,prevB2,prevC2 : real;
begin
  for l1:=0 to word1.number do begin
    for l2:=0 to word2.number do begin
      prevA1:= maxint;
      prevB1:= maxint;
      prevC1:= maxint;

      prevA2:=-maxint;
      prevB2:=-maxint;
      prevC2:=-maxint;
        
      if (l1>0)
        then begin
               prevA1:=matrix1[l1-1,l2  ]+weight(word1.letter[l1],empty);
               prevA2:=matrix2[l1-1,l2  ]+1;
             end;

      if (l1>0) and (l2>0)
        then begin
               prevB1:=matrix1[l1-1,l2-1]+weight(word1.letter[l1],word2.letter[l2]);
               prevB2:=matrix2[l1-1,l2-1]+1;
             end;

      if (l2>0)
        then begin
               prevC1:=matrix1[l1  ,l2-1]+weight(empty,word2.letter[l2]);
               prevC2:=matrix2[l1  ,l2-1]+1;
             end;

      check(matrix1[l1,l2],prevA1,prevB1,prevC1,
            matrix2[l1,l2],prevA2,prevB2,prevC2);

      if matrix1[l1,l2]= maxint then matrix1[l1,l2]:=0;
      if matrix2[l1,l2]=-maxint then matrix2[l1,l2]:=0;

      addtrace(prevA1,prevB1,prevC1,
               prevA2,prevB2,prevC2,l1,l2);
    end;
  end;

  distance0:=matrix1[word1.number,word2.number];
  
  proctrace(word1.number,word2.number,1,word1,word2,distance1,distance2);
end;

begin {main}
  getparameters;

  initcomp;
  initvc;
  initempty;

  rword(stringO, wordO);
  rword(stringY, wordY);
  
  change(wordO, wordO0);
  change(wordY, wordY0);
  
  determine(wordO0, wordY0, dist0, dist1, dist2);
  
  // write  (stdout, dist1      :0:6, ' ');
  // write  (stdout, dist2      :0:6, ' ');
  // writeln(stdout, dist1/dist2:0:6     );
end   {main}.
