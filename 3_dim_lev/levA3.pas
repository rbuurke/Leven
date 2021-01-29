{ (c) W. J. Heeringa 2011 }

{$R+ $B+}
program calculate_word_distance(input,output,stderr);
{Levenshtein/spectrograms}

const
  mm = 1000;
  mn = 1000;
  ml = 1000;
  mp = 1000;
  mt = 1000;

type
  dialecttype    = array[1..mn] of string[255];

  lexemetype     = array[1..3] of char;

  morphemetype   = array[1..3] of char;

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
                     lexeme    : lexemetype;
                     morpheme  : morphemetype;
                     letter    : array[1..80] of lettertype;
                     number    : integer;
                   end;

  listtype       = record
                     w_rd      : array[1..10] of wordtype;
                     number    : integer;
                   end;

  comptype       = record
                     number    : integer;
                     segment   : array[1..mp] of segmenttype;
                     matrix    : array[1..mp] of array[1..mp] of real;
                   end;

  segmenttype0   = record
                     sound     : integer;
                     long      : integer;
                     voice     : boolean;
                     apic      : boolean;
                     nasal     : boolean;
                     palat     : boolean;
                     velar     : boolean;
                   end;

  lettertype0    = record
                     segment   : array[1..2] of segmenttype0;
                   end;

  wordtype0      = record
                     lexeme    : integer;
                     morpheme  : integer;
                     letter    : array[1..80] of lettertype0;
                     number    : integer;
                   end;

  listtype0      = record
                     w_rd      : array[1..10] of wordtype0;
                     number    : integer;
                   end;

  matrixtype     = array[0..80,0..80,0..80] of real;

  tracetype      = array[0..80,0..80,0..80] of integer;

  alignrec       = record
                     word1     :  integer;
                     word2     :  integer;
                     word3     :  integer;
                     dist      :  real;
                   end;

  aligntype      = array[1..80] of alignrec;

  dtype          = array[1..mm] of real;

var
  m,n                          : integer;
  dialect,dialect0             : dialecttype;
  f                            : dtype;
  comp                         : comptype;
  vc                           : real;
  empty                        : lettertype0;
  matrix1,matrix2              : matrixtype;
  trace                        : tracetype;
  align                        : aligntype;
  fname,ffile,cfile,lfile,path : string[255];
  nr,id                        : integer;
  a_lex,a_mor,a_pho            : dtype;
  length,diphthong,part,format : integer;
  dnr,lnr                      : integer;
  fp_lex,fp_mor,fp_pho         : text;

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

procedure openw(var fp:text;name:string);
begin
  assign(fp,name);
  {$I-}
  rewrite(fp);
  {$I+}

  if IOResult <> 0
    then begin
           writeln(stderr,'Error opening file ',name);
           halt;
         end
    else {nothing};
end;

procedure usage;
begin
  writeln(stderr);
  writeln(stderr,'(c) W. J. Heeringa 2011');
  writeln(stderr);
  writeln(stderr,'Usage: ',ParamStr(0),' filfile frqfile cmpfile lblfile nr id 0|1|2 1|2 1-7 1|2');
  writeln(stderr);
  writeln(stderr,'filfile: file with files              ');
  writeln(stderr,'frqfile: file with frequencies        ');
  writeln(stderr,'cmpfile: file with segment comparisons');
  writeln(stderr,'lblfile: file with labels             ');
  writeln(stderr);
  writeln(stderr,'nr: number of standard languages');
  writeln(stderr,'id: index  of standard language ');
  writeln(stderr);
  writeln(stderr,'0: No   lengths');
  writeln(stderr,'1: Two  lengths');
  writeln(stderr,'2: Four lengths');
  writeln(stderr);
  writeln(stderr,'1: Diphthong is two segments');
  writeln(stderr,'2: Diphthong is one segment ');
  writeln(stderr);
  writeln(stderr,'1: Change                                        ');
  writeln(stderr,'2: Change due to convergence wrt standard language');
  writeln(stderr,'3: Change due to  neutrality wrt standard language');
  writeln(stderr,'4: Change due to  divergence wrt standard language');
  writeln(stderr);
  writeln(stderr,'5: Old distance to standard language');
  writeln(stderr,'6: New distance to standard language');
  writeln(stderr,'7: Convergence  to standard language');
  writeln(stderr);
  writeln(stderr,'1: links or vectors');
  writeln(stderr,'2: table           ');
  writeln(stderr);
  halt;
end;

procedure getparameters;
var
  code : integer;
begin
  if paramcount=10
    then begin
           fname:=paramstr(1);
           ffile:=paramstr(2);
           cfile:=paramstr(3);
           lfile:=paramstr(4);

           val(paramstr(5),nr,code);
           if code<>0
             then usage
             else {nothing};

           val(paramstr(6),id,code);
           if code<>0
             then usage
             else {nothing};

           if paramstr(7)='0' then length:=0 else
           if paramstr(7)='1' then length:=1 else
           if paramstr(7)='2' then length:=2 else usage;

           if paramstr(8)='1' then diphthong:=1 else
           if paramstr(8)='2' then diphthong:=2 else usage;

           if paramstr(9)='1' then part:=1 else
           if paramstr(9)='2' then part:=2 else
           if paramstr(9)='3' then part:=3 else
           if paramstr(9)='4' then part:=4 else
           if paramstr(9)='5' then part:=5 else
           if paramstr(9)='6' then part:=6 else
           if paramstr(9)='7' then part:=7 else usage;

           if paramstr(10)='1' then format:=1 else
           if paramstr(10)='2' then format:=2 else usage;
         end
    else usage;
end;

procedure initpath;
begin
  path:='../phon/';

  case length of
    0: begin
         path:=path+'leng0/';
       end;
    1: begin
         path:=path+'leng1/';
       end;
    2: begin
         path:=path+'leng2/';
       end;
  end;

  case diphthong of
    1: begin
         path:=path+'diph1/';
       end;
    2: begin
         path:=path+'diph2/';
       end;
  end;
end;

procedure initdialect;
var
  fp : text;
begin
  openr(fp,fname);

  n:=0;
  while not eof(fp) do begin
    inc(n);
    read(fp,dialect[n]);
    readln(fp);
  end;

  close(fp);
end;

procedure initdialect0;
var
  fp    : text;
  n     : integer;
  index : integer;
begin
  openr(fp,lfile);

  n:=0;
  while not eof(fp) do begin
    inc(n);
    readln(fp,index,dialect0[n]);
    while dialect0[n][1]=' ' do
      delete(dialect0[n],1,1);
  end;

  close(fp);
end;

function count(f:real):integer;
begin
  if f=0
    then count:=0
    else count:=1;
end;

procedure initfreq;
var
  fp : text;
  i  : integer;
begin
  openr(fp,ffile);

  m:=0;
  while not eof(fp) do begin
    read(fp,i);
    read(fp,f[i]);
    readln(fp);
    m:=m+count(f[i]);
  end;

  close(fp);
end;

function inlexemes(ch:char):boolean;
begin
  if ch='0'
    then inlexemes:=true
    else
  if ch='1'
    then inlexemes:=true
    else
  if ch='2'
    then inlexemes:=true
    else
  if ch='3'
    then inlexemes:=true
    else
  if ch='4'
    then inlexemes:=true
    else
  if ch='5'
    then inlexemes:=true
    else
  if ch='6'
    then inlexemes:=true
    else
  if ch='7'
    then inlexemes:=true
    else
  if ch='8'
    then inlexemes:=true
    else
  if ch='9'
    then inlexemes:=true
    else inlexemes:=false;
end;

function inmorphemes(ch:char):boolean;
begin
  if ch='0'
    then inmorphemes:=true
    else
  if ch='1'
    then inmorphemes:=true
    else
  if ch='2'
    then inmorphemes:=true
    else
  if ch='3'
    then inmorphemes:=true
    else
  if ch='4'
    then inmorphemes:=true
    else
  if ch='5'
    then inmorphemes:=true
    else
  if ch='6'
    then inmorphemes:=true
    else
  if ch='7'
    then inmorphemes:=true
    else
  if ch='8'
    then inmorphemes:=true
    else
  if ch='9'
    then inmorphemes:=true
    else inmorphemes:=false;
end;

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

procedure halt0(string0:string);
begin
  writeln(stderr,lnr:4,' ',string0);
  halt;
end;

procedure readchar(var fp:text;var ch:char);
begin
  if eoln(fp)
    then ch:='#'
    else read(fp,ch);
end;

procedure initlexeme(var lexeme:lexemetype);
var
  x : integer;
begin
  for x:=1 to 3 do begin
    lexeme[x]:=' ';
  end;
end;

procedure rlexeme(var fp:text;var lexeme:lexemetype;var ch:char);
begin
  initlexeme(lexeme);

  if inlexemes(ch)
    then begin
           lexeme[1]:=ch;
           readchar(fp,ch);
         end
    else halt0('numeral expected');

  if inlexemes(ch)
    then begin
           lexeme[2]:=ch;
           readchar(fp,ch);

           if inlexemes(ch)
             then begin
                    lexeme[3]:=ch;
                    readchar(fp,ch);
                  end
             else {nothing}
         end
    else {nothing};
end;

procedure initmorpheme(var morpheme:morphemetype);
var
  x : integer;
begin
  for x:=1 to 3 do begin
    morpheme[x]:=' ';
  end;
end;

procedure rmorpheme(var fp:text;var morpheme:morphemetype;var ch:char);
begin
  initmorpheme(morpheme);

  if inmorphemes(ch)
    then begin
           morpheme[1]:=ch;
           readchar(fp,ch);
         end
    else halt0('numeral expected');

  if inmorphemes(ch)
    then begin
           morpheme[2]:=ch;
           readchar(fp,ch);

           if inmorphemes(ch)
             then begin
                    morpheme[3]:=ch;
                    readchar(fp,ch);
                  end
             else {nothing}
         end
    else {nothing};
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

procedure rstress(var fp:text;var stress:stresstype;var ch:char);
begin
  if ch=''''
    then begin
           stress[1]:=ch;
           readchar(fp,ch);
         end
    else
  if ch='"'
    then begin
           stress[1]:=ch;
           readchar(fp,ch);
           if ch='"'
             then begin
                    stress[2]:=ch;
                    readchar(fp,ch);
                  end
             else {nothing}
         end
    else
  if ch='%'
    then begin
           stress[1]:=ch;
           readchar(fp,ch);
           if ch='%'
             then begin
                    stress[2]:=ch;
                    readchar(fp,ch);
                  end
             else {nothing}
         end
    else halt0(''' or % or " or ^ expected');
end;

procedure rhead(var fp:text;var head:headtype;var ch:char);
begin
  if inheads(ch)
    then begin
           head[1]:=ch;
           readchar(fp,ch);
           if ch='\'
             then begin
                    head[2]:=ch;
                    readchar(fp,ch);
                  end
             else {nothing};
           if ch='`'
             then begin
                    head[3]:=ch;
                    readchar(fp,ch);
                  end
             else {nothing}
         end
    else halt0('head expected')
end;

procedure rsupra(var fp:text;var supra:supratype;var ch:char);
begin
  if ch=':'
    then begin
           supra[1]:=ch;
           readchar(fp,ch);
         end
    else halt0(': expected');

  if ch='\'
    then begin
           supra[2]:=ch;
           readchar(fp,ch);
         end
    else {nothing};
end;

procedure rdiacritic(var fp:text;var diacritic:diacritictype;var ch:char);
begin
  if ch='_'
    then begin
           readchar(fp,ch);
           if indiacritics(ch)
             then begin
                    diacritic[1]:=ch;
                    readchar(fp,ch);
                  end
             else halt0('diacritic expected');

           if ch='\'
             then begin
                    diacritic[2]:=ch;
                    readchar(fp,ch);
                  end
             else {nothing};
         end
    else halt0('_ expected');
end;

procedure rsegment(var fp:text;var segment:segmenttype;var ch:char);
var
  p,d : integer;
begin
  p:=0;
  d:=0;

  rhead(fp,segment.head,ch);
  while insupras(ch) or (ch='_') do begin
    if insupras(ch)
      then begin
             inc(p);
             rsupra(fp,segment.supra[p],ch)
           end
      else begin
             inc(d);
             rdiacritic(fp,segment.diacritic[d],ch);
           end;
  end;
end;

procedure rletter(var fp:text;var letter:lettertype;var ch:char);
var
  s : integer;
begin
  s:=1;

  initletter(letter);
  if instresss(ch) or inheads(ch)
    then begin
           if instresss(ch)
             then rstress(fp,letter.stress,ch)
             else {nothing};
           rsegment(fp,letter.segment[s],ch)
         end
    else if ch='['
          then begin
                 readchar(fp,ch);
                 if instresss(ch)
                   then rstress(fp,letter.stress,ch)
                   else {nothing};
                 rsegment(fp,letter.segment[s],ch);

                 while inheads(ch) do begin
                   inc(s);
                   rsegment(fp,letter.segment[s],ch);
                 end;

                 if ch=']'
                   then readchar(fp,ch)
                   else halt0('] expected');
               end
          else halt0('stress or head or [ expected');
end;

procedure rword(var fp:text;var w_rd:wordtype;var ch:char);
var
  l : integer;
begin
  rlexeme(fp,w_rd.lexeme,ch);

  if ch=' '
    then readchar(fp,ch)
    else halt0('space expected');

  rmorpheme(fp,w_rd.morpheme,ch);

  if ch=' '
    then readchar(fp,ch)
    else halt0('space expected');

  l:=1;
  rletter(fp,w_rd.letter[l],ch);
  while instresss(ch) or inheads(ch) or (ch='[') do begin
    inc(l);
    rletter(fp,w_rd.letter[l],ch);
  end;
  w_rd.number:=l;
end;

procedure rlist(var fp:text; var list:listtype;var ch:char);
var
  w : integer;
begin
  w:=1;
  rword(fp,list.w_rd[w],ch);
  while ch=' ' do begin
    inc(w);
    readchar(fp,ch);
    if ch='/'
      then begin
             readchar(fp,ch);
             if ch=' '
               then begin
                      readchar(fp,ch);
                      rword(fp,list.w_rd[w],ch)
                    end
               else halt0('space expected')
           end
      else halt0('/ expected')
  end;
  list.number:=w;
end;

function lexeme(lexeme0:lexemetype):integer;
var
  value,value0,ec : integer;
begin
  val(lexeme0[1],value,ec);

  if lexeme0[2]<>' '
   then begin
          val(lexeme0[2],value0,ec);
          value:=(value*10) + value0;

          if lexeme0[3]<>' '
            then begin
                   val(lexeme0[3], value0,ec);
                   value:=(value*10) + value0;
                 end
            else {nothing}
         end
   else {nothing};

  lexeme:=value;
end;

function morpheme(morpheme0:morphemetype):integer;
var
  value,value0,ec : integer;
begin
  val(morpheme0[1],value,ec);

  if morpheme0[2]<>' '
   then begin
          val(morpheme0[2],value0,ec);
          value:=(value*10) + value0;

          if morpheme0[3]<>' '
            then begin
                   val(morpheme0[3], value0,ec);
                   value:=(value*10) + value0;
                 end
            else {nothing}
         end
   else {nothing};

  morpheme:=value;
end;

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

function long(segment:segmenttype):integer;
var
  found : boolean;
  d     : integer;
begin
  long:=2;

  d:=0;
  found:=false;

  while not found and (d<4) do begin
    d:=d+1;
    found:=(segment.diacritic[d]='X ');
  end;

  if found then long:=1;

  d:=0;
  found:=false;

  while not found and (d<4) do begin
    d:=d+1;
    found:=(segment.supra[d]=':\');
  end;

  if found then long:=3;

  d:=0;
  found:=false;

  while not found and (d<4) do begin
    d:=d+1;
    found:=(segment.supra[d]=': ');
  end;

  if found then long:=4;
end;

function voice(segment:segmenttype):boolean;
var
  found : boolean;
  d     : integer;
begin
  found:=false;

  d:=0;
  while not found and (d<4) do begin
    d:=d+1;
    found:=(segment.diacritic[d]='0 ') or
           (segment.diacritic[d]='v ');
  end;

  voice:=found;
end;

function apic(segment:segmenttype):boolean;
var
  found : boolean;
  d     : integer;
begin
  found:=false;

  d:=0;
  while not found and (d<4) do begin
    d:=d+1;
    found:=(segment.diacritic[d]='a ');
  end;

  apic:=found;
end;

function nasal(segment:segmenttype):boolean;
var
  found : boolean;
  d     : integer;
begin
  found:=false;

  d:=0;
  while not found and (d<4) do begin
    d:=d+1;
    found:=(segment.diacritic[d]='~ ');
  end;

  nasal:=found;
end;

function palat(segment:segmenttype):boolean;
var
  found : boolean;
  d     : integer;
begin
  found:=false;

  d:=0;
  while not found and (d<4) do begin
    d:=d+1;
    found:=(segment.diacritic[d]='j ');
  end;

  palat:=found;
end;

function velar(segment:segmenttype):boolean;
var
  found : boolean;
  d     : integer;
begin
  found:=false;

  d:=0;
  while not found and (d<4) do begin
    d:=d+1;
    found:=(segment.diacritic[d]='e ');
  end;

  velar:=found;
end;

procedure change(list:listtype;var list0:listtype0);
var
  w,l : integer;
begin
  for w:=1 to list.number do begin
    list0.w_rd[w].lexeme  :=lexeme  (list.w_rd[w].lexeme  );
    list0.w_rd[w].morpheme:=morpheme(list.w_rd[w].morpheme);

    for l:=1 to list.w_rd[w].number do begin
      list0.w_rd[w].letter[l].segment[1].sound:=sound(list.w_rd[w].letter[l].segment[1]);
      list0.w_rd[w].letter[l].segment[1].long :=long (list.w_rd[w].letter[l].segment[1]);
      list0.w_rd[w].letter[l].segment[1].voice:=voice(list.w_rd[w].letter[l].segment[1]);
      list0.w_rd[w].letter[l].segment[1].apic :=apic (list.w_rd[w].letter[l].segment[1]);
      list0.w_rd[w].letter[l].segment[1].nasal:=nasal(list.w_rd[w].letter[l].segment[1]);
      list0.w_rd[w].letter[l].segment[1].palat:=palat(list.w_rd[w].letter[l].segment[1]);
      list0.w_rd[w].letter[l].segment[1].velar:=velar(list.w_rd[w].letter[l].segment[1]);

      if list.w_rd[w].letter[l].segment[2].head='   '
        then begin
               list0.w_rd[w].letter[l].segment[2].sound:=list0.w_rd[w].letter[l].segment[1].sound;
               list0.w_rd[w].letter[l].segment[2].long :=list0.w_rd[w].letter[l].segment[1].long;
               list0.w_rd[w].letter[l].segment[2].sound:=list0.w_rd[w].letter[l].segment[1].sound;
               list0.w_rd[w].letter[l].segment[2].apic :=list0.w_rd[w].letter[l].segment[1].apic;
               list0.w_rd[w].letter[l].segment[2].nasal:=list0.w_rd[w].letter[l].segment[1].nasal;
               list0.w_rd[w].letter[l].segment[2].palat:=list0.w_rd[w].letter[l].segment[1].palat;
               list0.w_rd[w].letter[l].segment[2].velar:=list0.w_rd[w].letter[l].segment[1].velar;
             end
        else begin
               list0.w_rd[w].letter[l].segment[2].sound:=sound(list.w_rd[w].letter[l].segment[2]);
               list0.w_rd[w].letter[l].segment[2].long :=long (list.w_rd[w].letter[l].segment[2]);
               list0.w_rd[w].letter[l].segment[2].voice:=voice(list.w_rd[w].letter[l].segment[2]);
               list0.w_rd[w].letter[l].segment[2].apic :=apic (list.w_rd[w].letter[l].segment[2]);
               list0.w_rd[w].letter[l].segment[2].nasal:=nasal(list.w_rd[w].letter[l].segment[2]);
               list0.w_rd[w].letter[l].segment[2].palat:=palat(list.w_rd[w].letter[l].segment[2]);
               list0.w_rd[w].letter[l].segment[2].velar:=velar(list.w_rd[w].letter[l].segment[2]);
             end;
    end;

    list0.w_rd[w].number:=list.w_rd[w].number;
  end;

  list0.number:=list.number;
end;

procedure rline(var fp:text;var list0:listtype0);
var
  list : listtype;
  ch   : char;
begin
  readchar(fp,ch);

  if inlexemes(ch)
    then rlist(fp,list,ch)
    else if ch='|'
           then begin
                  list.number:=0;
                  readchar(fp,ch);
                end
           else halt0('numeral or | expected');

  if ch='#'
    then readln(fp)
    else halt0('eoln expected');

  change(list,list0);
end;

procedure initcomp;
var
  fp  : text;
  i,j : integer;
  ch  : char;
  v   : real;
begin
  openr(fp,cfile);

  readln(fp,comp.number);

  for i:=1 to comp.number do begin
    initsegment(comp.segment[i]);
    read(fp,ch);
    rsegment(fp,comp.segment[i],ch);
    readln(fp);
    comp.matrix[i,i]:=0;
  end;

  for i:=2 to comp.number do begin
    for j:=1 to (i-1) do begin
      readln(fp,v);
      comp.matrix[i,j]:=ln(v+1);
      comp.matrix[j,i]:=ln(v+1);
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

procedure initempty;
var
  s : integer;
begin
  for s:=1 to 2 do begin
    empty.segment[s].sound:=1;
    empty.segment[s].long :=0;
    empty.segment[s].voice:=false;
    empty.segment[s].apic :=false;
    empty.segment[s].nasal:=false;
    empty.segment[s].palat:=false;
    empty.segment[s].velar:=false;
  end;
end;

function vowel1(sound:integer):boolean;
begin
  vowel1:=
     ((sound>=02) and (sound<=29)) or
      (sound =82) or  (sound =88);
end;

function vowel10(sound:integer):boolean;
begin
  vowel10:=
     ((sound>=02) and (sound<=29))
end;

function consonant1(sound:integer):boolean;
begin
  consonant1:=
     ((sound>=30) and (sound<=88)) or
      (sound =02) or  (sound =07);
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

function voiced(var sound:integer):boolean;
begin
  voiced:=true;

  case sound of
    30 : sound:=31;
    32 : sound:=33;
    34 : sound:=35;
    36 : sound:=37;
    38 : sound:=39;
    40 : sound:=41;
    55 : sound:=56;
    57 : sound:=58;
    59 : sound:=60;
    61 : sound:=62;
    63 : sound:=64;
    65 : sound:=66;
    67 : sound:=68;
    69 : sound:=70;
    71 : sound:=72;
    73 : sound:=74;
    75 : sound:=76;

    else voiced:=false;
  end;
end;

function voiceless(var sound:integer):boolean;
begin
  voiceless:=true;

  case sound of
    31 : sound:=30;
    33 : sound:=32;
    35 : sound:=34;
    37 : sound:=36;
    39 : sound:=38;
    41 : sound:=40;
    56 : sound:=55;
    58 : sound:=57;
    60 : sound:=59;
    62 : sound:=61;
    64 : sound:=63;
    66 : sound:=65;
    68 : sound:=67;
    70 : sound:=69;
    72 : sound:=71;
    74 : sound:=73;
    76 : sound:=75;

    else voiceless:=false;
  end;
end;

function opposite(sound:integer):integer;
begin
  if voiced   (sound) then opposite:=sound else
  if voiceless(sound) then opposite:=sound else opposite:=sound;
end;

function weight00(segment1,segment2:segmenttype0):real;
var
  d : real;
  k : integer;
begin
  if (    vowel1(segment1.sound) and     vowel1(segment2.sound)) or
     (consonant1(segment1.sound) and consonant1(segment2.sound)) or

     (    vowel2(segment1.sound) and consonant2(segment2.sound)) or
     (consonant2(segment1.sound) and     vowel2(segment2.sound)) or

     (    vowel1(segment1.sound) and     empty0(segment2.sound)) or
     (    empty0(segment1.sound) and     vowel1(segment2.sound)) or
     (consonant1(segment1.sound) and     empty0(segment2.sound)) or
     (    empty0(segment1.sound) and consonant1(segment2.sound)) or
     (    empty0(segment1.sound) and     empty0(segment2.sound))

    then begin
           d:=comp.matrix[segment1.sound,segment2.sound];
           k:=1;

           if (segment1.voice) and (not segment2.voice)
             then begin
                    d:=d+comp.matrix[opposite(segment1.sound),segment2.sound];
                    k:=k+1;
                  end
             else
           if (not segment1.voice) and (segment2.voice)
             then begin
                    d:=d+comp.matrix[segment1.sound,opposite(segment2.sound)];
                    k:=k+1;
                  end
             else
           if (segment1.voice) and (segment2.voice)
             then begin
                    d:=d+comp.matrix[opposite(segment1.sound),opposite(segment2.sound)];
                    k:=k+1;
                  end
             else {nothing};

           if (segment1.apic) and (not segment2.apic)
             then if (segment1.sound=61)
                    then begin
                           d:=d+comp.matrix[65,segment2.sound];
                           k:=k+1;
                         end
                    else
                  if (segment1.sound=62)
                    then begin
                           d:=d+comp.matrix[66,segment2.sound];
                           k:=k+1;
                         end
                    else halt0('apical applied to wrong sound')
             else
           if (not segment1.apic) and (segment2.apic)
             then if (segment2.sound=61)
                    then begin
                           d:=d+comp.matrix[segment1.sound,65];
                           k:=k+1;
                         end
                    else
                  if (segment2.sound=62)
                    then begin
                           d:=d+comp.matrix[segment1.sound,66];
                           k:=k+1;
                         end
                    else halt0('apical applied to wrong sound')
             else
           if (segment1.apic) and (segment2.apic)
             then if (segment1.sound=61) and (segment2.sound=62)
                    then begin
                           d:=d+comp.matrix[65,66];
                           k:=k+1;
                         end
                    else
                  if (segment1.sound=62) and (segment2.sound=61)
                    then begin
                           d:=d+comp.matrix[66,65];
                           k:=k+1;
                         end
                    else
                  if (segment1.sound=61) and (segment2.sound=61)
                    then begin
                           d:=d+comp.matrix[65,65];
                           k:=k+1;
                         end
                    else
                  if (segment1.sound=62) and (segment2.sound=62)
                    then begin
                           d:=d+comp.matrix[66,66];
                           k:=k+1;
                         end
                    else halt0('apical applied to wrong sound')
             else {nothing};

           if (segment1.nasal) and (not segment2.nasal)
             then begin
                    d:=d+comp.matrix[45,segment2.sound];
                    k:=k+1;
                  end
             else
           if (not segment1.nasal) and (segment2.nasal)
             then begin
                    d:=d+comp.matrix[segment1.sound,45];
                    k:=k+1;
                  end
             else
           if (segment1.nasal) and (segment2.nasal)
             then begin
                    d:=d+(comp.matrix[45,45]);
                    k:=k+1;
                  end
             else {nothing};

           if (segment1.palat) and (not segment2.palat)
             then begin
                    d:=d+comp.matrix[82,segment2.sound];
                    k:=k+1;
                  end
             else
           if (not segment1.palat) and (segment2.palat)
             then begin
                    d:=d+comp.matrix[segment1.sound,82];
                    k:=k+1;
                  end
             else
           if (segment1.palat) and (segment2.palat)
             then begin
                    d:=d+(comp.matrix[82,82]);
                    k:=k+1;
                  end
             else {nothing};

           if (segment1.velar) and (not segment2.velar)
             then begin
                    d:=d+comp.matrix[07,segment2.sound];
                    k:=k+1;
                  end
             else
           if (not segment1.velar) and (segment2.velar)
             then begin
                    d:=d+comp.matrix[segment1.sound,07];
                    k:=k+1;
                  end
             else
           if (segment1.velar) and (segment2.velar)
             then begin
                    d:=d+(comp.matrix[07,07]);
                    k:=k+1;
                  end
             else {nothing};

           weight00:=d/k;
         end
    else weight00:=maxint;
end;

function weight0(letter1,letter2:lettertype0):real;
begin
  weight0:=(weight00(letter1.segment[1],letter2.segment[1])+
            weight00(letter1.segment[2],letter2.segment[2]))/2;
end;

function weight(letter1,letter2,letter3:lettertype0):real;
begin
  weight:=weight0(letter1,letter2)+
          weight0(letter1,letter3)+
          weight0(letter2,letter3);
end;

function min(a,b,c,d,e,f,g:real):real;
begin
  if (a<b) and (a<c) and (a<d) and (a<e) and (a<f) and (a<g)
    then min:=a
    else

  if (b<c) and (b<d) and (b<e) and (b<f) and (b<g)
    then min:=b
    else

  if (c<d) and (c<e) and (c<f) and (c<g)
    then min:=c
    else

  if (d<e) and (d<f) and (d<g)
    then min:=d
    else

  if (e<f) and (e<g)
    then min:=e
    else

  if (f<g)
    then min:=f
    else min:=g;
end;

function max(a,b,c,d,e,f,g:real):real;
begin
  if (a>b) and (a>c) and (a>d) and (a>e) and (a>f) and (a>g)
    then max:=a
    else

  if (b>c) and (b>d) and (b>e) and (b>f) and (b>g)
    then max:=b
    else

  if (c>d) and (c>e) and (c>f) and (c>g)
    then max:=c
    else

  if (d>e) and (d>f) and (d>g)
    then max:=d
    else

  if (e>f) and (e>g)
    then max:=e
    else

  if (f>g)
    then max:=f
    else max:=g;
end;

procedure check(var value1:real;prevA1,prevB1,prevC1,prevD1,prevE1,prevF1,prevG1:real;
                var value2:real;prevA2,prevB2,prevC2,prevD2,prevE2,prevF2,prevG2:real);
begin
  value1:=min(prevA1,prevB1,prevC1,prevD1,prevE1,prevF1,prevG1);

  if prevA1<>value1
    then prevA2:=-maxint;

  if prevB1<>value1
    then prevB2:=-maxint;

  if prevC1<>value1
    then prevC2:=-maxint;

  if prevD1<>value1
    then prevD2:=-maxint;

  if prevE1<>value1
    then prevE2:=-maxint;

  if prevF1<>value1
    then prevF2:=-maxint;

  if prevG1<>value1
    then prevG2:=-maxint;

  value2:=max(prevA2,prevB2,prevC2,prevD2,prevE2,prevF2,prevG2);
end;

procedure addtrace(prevA1,prevB1,prevC1,prevD1,prevE1,prevF1,prevG1,
                   prevA2,prevB2,prevC2,prevD2,prevE2,prevF2,prevG2:real;l1,l2,l3:integer);
var
  pointer : integer;
begin
  pointer:=0;

  if (prevA1=matrix1[l1,l2,l3]) and (prevA2=matrix2[l1,l2,l3])
    then pointer:=pointer+128;

  if (prevB1=matrix1[l1,l2,l3]) and (prevB2=matrix2[l1,l2,l3])
    then pointer:=pointer+64;

  if (prevC1=matrix1[l1,l2,l3]) and (prevC2=matrix2[l1,l2,l3])
    then pointer:=pointer+32;

  if (prevD1=matrix1[l1,l2,l3]) and (prevD2=matrix2[l1,l2,l3])
    then pointer:=pointer+16;

  if (prevE1=matrix1[l1,l2,l3]) and (prevE2=matrix2[l1,l2,l3])
    then pointer:=pointer+8;

  if (prevF1=matrix1[l1,l2,l3]) and (prevF2=matrix2[l1,l2,l3])
    then pointer:=pointer+4;

  if (prevG1=matrix1[l1,l2,l3]) and (prevG2=matrix2[l1,l2,l3])
    then pointer:=pointer+2;

  trace[l1,l2,l3]:=pointer;
end;

procedure addalign(l1,l2,l3:integer;d:real;a:integer);
begin
  align[a].word1:=l1;
  align[a].word2:=l2;
  align[a].word3:=l3;
  align[a].dist :=d ;
end;

procedure procalign(a:integer;word1,word2,word3:wordtype0;var distance1,distance2:real);
var
  l                 : integer;
  disXX,oldXS,newXS : real;
begin
  distance1:=0;
  distance2:=0;

  for l:=(a-1) downto 1 do begin
    { distance between oldX and newX }

    if (align[l].word1<>0) and (align[l].word2<>0)
      then disXX:=weight0(word1.letter[align[l].word1],word2.letter[align[l].word2])
      else

    if (align[l].word1<>0) and (align[l].word2= 0)
      then disXX:=weight0(word1.letter[align[l].word1],empty)
      else

    if (align[l].word1= 0) and (align[l].word2<>0)
      then disXX:=weight0(empty,word2.letter[align[l].word2])
      else disXX:=0;

    { distance between oldX and S }

    if (align[l].word1<>0) and (align[l].word3<>0)
      then oldXS:=weight0(word1.letter[align[l].word1],word3.letter[align[l].word3])
      else

    if (align[l].word1<>0) and (align[l].word3= 0)
      then oldXS:=weight0(word1.letter[align[l].word1],empty)
      else

    if (align[l].word1= 0) and (align[l].word3<>0)
      then oldXS:=weight0(empty,word3.letter[align[l].word3])
      else oldXS:=0;

    { distance between newX and S }

    if (align[l].word2<>0) and (align[l].word3<>0)
      then newXS:=weight0(word2.letter[align[l].word2],word3.letter[align[l].word3])
      else

    if (align[l].word2<>0) and (align[l].word3= 0)
      then newXS:=weight0(word2.letter[align[l].word2],empty)
      else

    if (align[l].word2= 0) and (align[l].word3<>0)
      then newXS:=weight0(empty,word3.letter[align[l].word3])
      else newXS:=0;


    if  part=1
      then distance1:=distance1+disXX
      else

    if (part=2) and (oldXS>newXS)
      then distance1:=distance1+disXX
      else

    if (part=3) and (oldXS=newXS)
      then distance1:=distance1+disXX
      else

    if (part=4) and (oldXS<newXS)
      then distance1:=distance1+disXX
      else

    if  part=5
      then distance1:=distance1+oldXS
      else

    if  part=6
      then distance1:=distance1+newXS
      else

    if  part=7
      then distance1:=distance1+(oldXS-newXS)
      else {nothing};

    distance2:=distance2+1;
  end;

  distance1:=(distance1*100)/vc;
end;

procedure proctrace(l1,l2,l3,a:integer;word1,word2,word3:wordtype0;var distance1,distance2:real);
var
  pointer : integer;
begin
  if (l1>0) or (l2>0) or (l3>0)
    then begin
           pointer:=trace[l1,l2,l3];

           if (pointer>=128) and (l1>0)
             then begin
                    addalign(l1,0,0,matrix1[l1,l2,l3]-matrix1[l1-1,l2,l3],a);
                    proctrace(l1-1,l2,l3,a+1,word1,word2,word3,distance1,distance2);
                    pointer:=pointer-128;
                  end;

           if (pointer>=64) and (l2>0)
             then begin
                    addalign(0,l2,0,matrix1[l1,l2,l3]-matrix1[l1,l2-1,l3],a);
                    proctrace(l1,l2-1,l3,a+1,word1,word2,word3,distance1,distance2);
                    pointer:=pointer-64;
                  end;

           if (pointer>=32) and (l3>0)
             then begin
                    addalign(0,0,l3,matrix1[l1,l2,l3]-matrix1[l1,l2,l3-1],a);
                    proctrace(l1,l2,l3-1,a+1,word1,word2,word3,distance1,distance2);
                    pointer:=pointer-32;
                  end;

           if (pointer>=16) and (l1>0) and (l2>0)
             then begin
                    addalign(l1,l2,0,matrix1[l1,l2,l3]-matrix1[l1-1,l2-1,l3],a);
                    proctrace(l1-1,l2-1,l3,a+1,word1,word2,word3,distance1,distance2);
                    pointer:=pointer-16;
                  end;

           if (pointer>=8) and (l1>0) and (l3>0)
             then begin
                    addalign(l1,0,l3,matrix1[l1,l2,l3]-matrix1[l1-1,l2,l3-1],a);
                    proctrace(l1-1,l2,l3-1,a+1,word1,word2,word3,distance1,distance2);
                    pointer:=pointer-8;
                  end;

           if (pointer>=4) and (l2>0) and (l3>0)
             then begin
                    addalign(0,l2,l3,matrix1[l1,l2,l3]-matrix1[l1,l2-1,l3-1],a);
                    proctrace(l1,l2-1,l3-1,a+1,word1,word2,word3,distance1,distance2);
                    pointer:=pointer-4;
                  end;

           if (pointer>=2) and (l1>0) and (l2>0) and (l3>0)
             then begin
                    addalign(l1,l2,l3,matrix1[l1,l2,l3]-matrix1[l1-1,l2-1,l3-1],a);
                    proctrace(l1-1,l2-1,l3-1,a+1,word1,word2,word3,distance1,distance2);
                    pointer:=pointer-2;
                  end;
         end
    else procalign(a,word1,word2,word3,distance1,distance2);
end;

procedure determine(word1,word2,word3:wordtype0;var distance0,distance1,distance2:real);
var
  l1,l2,l3                                         : integer;
  prevA1,prevB1,prevC1,prevD1,prevE1,prevF1,prevG1 : real;
  prevA2,prevB2,prevC2,prevD2,prevE2,prevF2,prevG2 : real;
begin
  if (word1.lexeme=word2.lexeme) and (word1.morpheme=word2.morpheme) and
     (word1.lexeme=word3.lexeme) and (word1.morpheme=word3.morpheme) and
     (word2.lexeme=word3.lexeme) and (word2.morpheme=word3.morpheme)
    then begin
           for l1:=0 to word1.number do begin
             for l2:=0 to word2.number do begin
               for l3:=0 to word3.number do begin
                 prevA1:= maxint;
                 prevB1:= maxint;
                 prevC1:= maxint;
                 prevD1:= maxint;
                 prevE1:= maxint;
                 prevF1:= maxint;
                 prevG1:= maxint;

                 prevA2:=-maxint;
                 prevB2:=-maxint;
                 prevC2:=-maxint;
                 prevD2:=-maxint;
                 prevE2:=-maxint;
                 prevF2:=-maxint;
                 prevG2:=-maxint;

                 if (l1>0)
                   then begin
                          prevA1:=matrix1[l1-1,l2,l3]+weight(word1.letter[l1],empty,empty);
                          prevA2:=matrix2[l1-1,l2,l3]+1;
                        end;

                 if (l2>0)
                   then begin
                          prevB1:=matrix1[l1,l2-1,l3]+weight(empty,word2.letter[l2],empty);
                          prevB2:=matrix2[l1,l2-1,l3]+1;
                        end;

                 if (l3>0)
                   then begin
                          prevC1:=matrix1[l1,l2,l3-1]+weight(empty,empty,word3.letter[l3]);
                          prevC2:=matrix2[l1,l2,l3-1]+1;
                        end;

                 if (l1>0) and (l2>0)
                   then begin
                          prevD1:=matrix1[l1-1,l2-1,l3]+weight(word1.letter[l1],word2.letter[l2],empty);
                          prevD2:=matrix2[l1-1,l2-1,l3]+1;
                        end;

                 if (l1>0) and (l3>0)
                   then begin
                          prevE1:=matrix1[l1-1,l2,l3-1]+weight(word1.letter[l1],empty,word3.letter[l3]);
                          prevE2:=matrix2[l1-1,l2,l3-1]+1;
                        end;

                 if (l2>0) and (l3>0)
                   then begin
                          prevF1:=matrix1[l1,l2-1,l3-1]+weight(empty,word2.letter[l2],word3.letter[l3]);
                          prevF2:=matrix2[l1,l2-1,l3-1]+1;
                        end;

                 if (l1>0) and (l2>0) and (l3>0)
                   then begin
                          prevG1:=matrix1[l1-1,l2-1,l3-1]+weight(word1.letter[l1],word2.letter[l2],word3.letter[l3]);
                          prevG2:=matrix2[l1-1,l2-1,l3-1]+1;
                        end;

                 check(matrix1[l1,l2,l3],prevA1,prevB1,prevC1,prevD1,prevE1,prevF1,prevG1,
                       matrix2[l1,l2,l3],prevA2,prevB2,prevC2,prevD2,prevE2,prevF2,prevG2);

                 if matrix1[l1,l2,l3]= maxint then matrix1[l1,l2,l3]:=0;
                 if matrix2[l1,l2,l3]=-maxint then matrix2[l1,l2,l3]:=0;

                 addtrace(prevA1,prevB1,prevC1,prevD1,prevE1,prevF1,prevG1,
                          prevA2,prevB2,prevC2,prevD2,prevE2,prevF2,prevG2,l1,l2,l3);
               end;
             end;
           end;

           distance0:=matrix1[word1.number,word2.number,word3.number];
           proctrace(word1.number,word2.number,word3.number,1,word1,word2,word3,distance1,distance2);
         end
    else distance0:=1.7e+38;
end;

procedure calculate(i,j,k:integer;var a_lex,a_mor,a_pho:real);
var
  list1,list2,list3          : listtype0;
  dialect1,dialect2,dialect3 : text;
  w1,w2,w3                   : integer;
  matrix0,matrix1,matrix2    : array[1..10] of array[1..10] of array[1..10] of real;
  c1,c2,c3                   : array[0..10] of integer;
  h1,h2,h3                   : integer;
  d00                        : real;
  w10,w20,w30                : integer;
  d_pho,d_len,a_len          : real;
  n_pho,n_len                : integer;
begin
  openr(dialect1,path+dialect[i]);
  openr(dialect2,path+dialect[j]);
  openr(dialect3,path+dialect[k]);

  lnr:=0; a_pho:=0; a_len:=0;
  while (not eof(dialect1)) and (not eof(dialect2)) and (not eof(dialect3)) do begin
    lnr:=lnr+1;

    rline(dialect1,list1);
    rline(dialect2,list2);
    rline(dialect3,list3);

    if (f[lnr]<>0)
      then if (list1.number>0) and (list2.number>0) and (list3.number>0)
             then begin
                    for w1:=1 to list1.number do
                      for w2:=1 to list2.number do
                        for w3:=1 to list3.number do
                          determine(list1.w_rd[w1],list2.w_rd[w2],list3.w_rd[w3],matrix0[w1,w2,w3],matrix1[w1,w2,w3],matrix2[w1,w2,w3]);

                    for w1:=1 to list1.number do
                      c1[w1]:=list2.number*list3.number;
                    for w2:=1 to list2.number do
                      c2[w2]:=list1.number*list3.number;
                    for w3:=1 to list3.number do
                      c3[w3]:=list1.number*list2.number;

                    d_pho:=0;n_pho:=0;
                    d_len:=0;n_len:=0;

                    for h1:=1 to list1.number do begin
                      for h2:=1 to list2.number do begin
                        for h3:=1 to list3.number do begin

                          d00:=1.7e+38;
                          for w1:=1 to list1.number do begin
                            for w2:=1 to list2.number do begin
                              for w3:=1 to list3.number do begin
                                if (matrix0[w1,w2,w3]<d00) and (c1[w1]>0) and (c2[w2]>0) and (c3[w3]>0)
                                  then begin
                                         d00:=matrix0[w1,w2,w3];
                                         w10:=w1;
                                         w20:=w2;
                                         w30:=w3;
                                       end
                                  else {nothing};
                              end;
                            end;
                          end;

                          if (d00<1.7e+38)
                            then begin
                                   d_pho:=d_pho+matrix1[w10,w20,w30];
                                   d_len:=d_len+matrix2[w10,w20,w30];

                                   n_pho:=n_pho+1;
                                   n_len:=n_len+1;
                                 end
                            else {nothing};

                          c1[w10]:=c1[w10]-1;
                          c2[w20]:=c2[w20]-1;
                          c3[w30]:=c3[w30]-1;
                        end;
                      end;
                    end;

                    if (n_pho>0)
                      then begin
                             d_pho:=(d_pho/n_pho)*f[lnr];
                             d_len:=(d_len/n_len)*f[lnr];

                             a_pho:=a_pho+d_pho;
                             a_len:=a_len+d_len;
                           end
                      else {nothing};
                  end
             else {nothing}
      else {nothing};
  end;

  close(dialect1);
  close(dialect2);
  close(dialect3);

  a_pho:=a_pho/a_len;
end;

procedure process1(var fp:text;a:dtype);
var
  dnr                       : integer;
  smallest,largest,distance : real;
begin
  smallest:=1.7E+38;
  largest :=0;

  for dnr:=(nr+1) to (nr+((n-nr) div 2)) do begin
    if a[dnr]<smallest
      then smallest:=a[dnr]
      else {nothing};

    if a[dnr]>largest
      then largest:=a[dnr]
      else {nothing};
  end;

  writeln(fp,3);

  for dnr:=(nr+1) to (nr+((n-nr) div 2)) do begin
    writeln(fp,dialect0[dnr]);

    if (smallest<largest)
      then distance:=(a[dnr]-smallest)/(largest-smallest)
      else distance:=0;

    writeln(fp,1-distance);
    writeln(fp,1-distance);
    writeln(fp,1);
  end;
end;

procedure process2(var fp:text;a:dtype);
var
  dnr     : integer;
  largest : real;
begin
  largest:=0;
  for dnr:=(nr+1) to (nr+((n-nr) div 2)) do begin
    if abs(a[dnr])>largest
      then largest:=abs(a[dnr])
      else {nothing};
  end;

  writeln(fp,3);
  for dnr:=(nr+1) to (nr+((n-nr) div 2)) do begin
    writeln(fp,dialect0[dnr]);

    if (a[dnr]>0)
      then begin
             writeln(fp,1);
             writeln(fp,1-abs(a[dnr]/largest));
             writeln(fp,1-abs(a[dnr]/largest));
           end
      else

     if (a[dnr]=0)
      then begin
             writeln(fp,1);
             writeln(fp,1);
             writeln(fp,1);
           end
      else

    if (a[dnr]<0)
      then begin
             writeln(fp,1-abs(a[dnr]/largest));
             writeln(fp,1-abs(a[dnr]/largest));
             writeln(fp,1);
           end
      else {nothing};
  end;
end;

procedure process3(var fp:text;a:dtype);
var
  dnr : integer;
begin
  for dnr:=(nr+1) to (nr+((n-nr) div 2)) do begin
    writeln(fp,a[dnr]);
  end;
end;

begin {main}
  getparameters;

  initpath;
  initdialect;
  initdialect0;
  initfreq;

  initcomp;
  initvc;
  initempty;

  writeln(stderr);
  for dnr:=(nr+1) to (nr+((n-nr) div 2)) do begin
    writeln(stderr,'Processing ',dialect[dnr],' and ',dialect[dnr + ((n-nr) div 2)],' and ',dialect[id]);
    calculate(dnr,dnr + ((n-nr) div 2),id,a_lex[dnr],a_mor[dnr],a_pho[dnr]);
  end;
  writeln(stderr);

  if format=1
    then begin
          {openw(fp_lex,'data.lex.ve3');}
          {openw(fp_mor,'data.mor.ve3');}
           openw(fp_pho,'data.pho.ve3');

           case part of
             1 : process1(fp_pho,a_pho);
             2 : process1(fp_pho,a_pho);
             3 : process1(fp_pho,a_pho);
             4 : process1(fp_pho,a_pho);
             5 : process1(fp_pho,a_pho);
             6 : process1(fp_pho,a_pho);
             7 : process2(fp_pho,a_pho);
           end
         end
    else

  if format=2
    then begin
          {openw(fp_lex,'table_lex.txt');}
          {openw(fp_mor,'table_mor.txt');}
           openw(fp_pho,'table_pho.txt');

           case part of
             1 : process3(fp_pho,a_pho);
             2 : process3(fp_pho,a_pho);
             3 : process3(fp_pho,a_pho);
             4 : process3(fp_pho,a_pho);
             5 : process3(fp_pho,a_pho);
             6 : process3(fp_pho,a_pho);
             7 : process3(fp_pho,a_pho);
           end
         end
    else {nothing};

 {close(fp_lex);}
 {close(fp_mor);}
  close(fp_pho);
end   {main}.
