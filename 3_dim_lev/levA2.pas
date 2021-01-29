{ (c) W. J. Heeringa 2008 }

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

  comptype1      = record
                     number    : integer;
                     segment   : array[1..mp] of segmenttype;
                     matrix    : array[1..mp] of array[1..mp] of real;
                   end;

  comptype2      = record
                     number    : integer;
                     stress    : array[1..mt] of stresstype;
                     matrix    : array[1..mt] of array[1..mt] of real;
                   end;

  segmenttype0   = record
                     sound     : integer;
                     voice     : boolean;
                     apic      : boolean;
                     nasal     : boolean;
                     palat     : boolean;
                     velar     : boolean;
                   end;

  lettertype0    = record
                     stress    : integer;
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

  wltype         = array[1..mm] of array[0..ml] of real;

  matrixtype     = array[0..80,0..80] of real;

  tracetype      = array[0..80,0..80] of integer;

  alignrec       = record
                     word1     :  integer;
                     word2     :  integer;
                     dist      :  real;
                   end;

  aligntype      = array[1..80] of alignrec;

  dtype          = array[1..mm] of real;

var
  m,n                                    : integer;
  i,j                                    : integer;
  dialect                                : dialecttype;
  f                                      : dtype;
  m0_lex,m0_mor,m0_pho,m0_ton,m0_len,m00 : real;
  comp1                                  : comptype1;
  comp2                                  : comptype2;
  vc1,vc2                                : real;
  matrix1,matrix2,matrix3                : matrixtype;
  trace                                  : tracetype;
  align                                  : aligntype;
  wl1,wl2                                : wltype;
  d_lex,d_mor,d_pho,d_ton,d_len          : dtype;
  fname,ffile,cfile1,cfile2,path         : string[255];
  length,diphthong,giw1,giw2,part        : integer;
  lnr                                    : integer;
  fp_lex,fp_mor,fp_pho,fp_ton            : text;

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
  writeln(stderr,'(c) W. J. Heeringa 2008');
  writeln(stderr);
  writeln(stderr,'Usage: ',ParamStr(0),' filfile frqfile cmpfile1 cmpfile2 0|1|2 1|2 1|2 1|2 [1-7]');
  writeln(stderr);
  writeln(stderr,'filfile : file with files              ');
  writeln(stderr,'frqfile : file with frequencies        ');
  writeln(stderr,'cmpfile1: file with segment comparisons');
  writeln(stderr,'cmpfile2: file with stress  comparisons');
  writeln(stderr);
  writeln(stderr,'0: No   lengths');
  writeln(stderr,'1: Two  lengths');
  writeln(stderr,'2: Four lengths');
  writeln(stderr);
  writeln(stderr,'1: Diphthong is two segments');
  writeln(stderr,'2: Diphthong is one segment ');
  writeln(stderr);
  writeln(stderr,'1: Lexis      excl. GIW');
  writeln(stderr,'2: Lexis      incl. GIW');
  writeln(stderr);
  writeln(stderr,'1: Morphology excl. GIW');
  writeln(stderr,'2: Morphology incl. GIW');
  writeln(stderr);
  writeln(stderr,'1: Only vowels                               ');
  writeln(stderr,'2: Only vowel     vs. vowel     substitutions');
  writeln(stderr,'3: Only vowel     indels                     ');
  writeln(stderr,'4: Only consonants                           ');
  writeln(stderr,'5: Only consonant vs. consonant substitutions');
  writeln(stderr,'6: Only consonant indels                     ');
  writeln(stderr,'7: Only schwa     vs. sonorant  substitutions');
  writeln(stderr);
  halt;
end;

procedure getparameters;
begin
  if (paramcount=8) or (paramcount=9)
    then begin
           fname :=paramstr(1);
           ffile :=paramstr(2);
           cfile1:=paramstr(3);
           cfile2:=paramstr(4);

           if paramstr(5)='0' then length:=0 else
           if paramstr(5)='1' then length:=1 else
           if paramstr(5)='2' then length:=2 else usage;

           if paramstr(6)='1' then diphthong:=1 else
           if paramstr(6)='2' then diphthong:=2 else usage;

           if paramstr(7)='1' then giw1:=1 else
           if paramstr(7)='2' then giw1:=2 else usage;

           if paramstr(8)='1' then giw2:=1 else
           if paramstr(8)='2' then giw2:=2 else usage;
         end
    else usage;

  if (paramcount=8)
    then part:=0
    else begin
           if paramstr(9)='1' then part:=1 else
           if paramstr(9)='2' then part:=2 else
           if paramstr(9)='3' then part:=3 else
           if paramstr(9)='4' then part:=4 else
           if paramstr(9)='5' then part:=5 else
           if paramstr(9)='6' then part:=6 else
           if paramstr(9)='7' then part:=7 else usage;
         end;
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
  s,p,d : integer;
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

function tone(letter:lettertype):integer;
var
  t     : integer;
  found : boolean;
begin
  t:=0;
  found:=false;

  while (not found) and (t<comp2.number) do begin
    t:=t+1;
    found:=(letter.stress=comp2.stress[t]);
  end;

  if not found
    then begin
           writeln(stderr,letter.stress,' not found in comparison table');
           halt;
         end
    else {nothing};

  tone:=t;
end;

function sound(segment:segmenttype):integer;
var
  t     : integer;
  found : boolean;
begin
  t:=0;
  found:=false;

  while (not found) and (t<comp1.number) do begin
    t:=t+1;
    found:=(segment.head=comp1.segment[t].head);
  end;

  if not found
    then begin
           writeln(stderr,segment.head,' not found in comparison table');
           halt;
         end
    else {nothing};

  sound:=t;
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
      list0.w_rd[w].letter[l].stress:=tone(list.w_rd[w].letter[l]);

      list0.w_rd[w].letter[l].segment[1].sound:=sound(list.w_rd[w].letter[l].segment[1]);
      list0.w_rd[w].letter[l].segment[1].voice:=voice(list.w_rd[w].letter[l].segment[1]);
      list0.w_rd[w].letter[l].segment[1].apic :=apic (list.w_rd[w].letter[l].segment[1]);
      list0.w_rd[w].letter[l].segment[1].nasal:=nasal(list.w_rd[w].letter[l].segment[1]);
      list0.w_rd[w].letter[l].segment[1].palat:=palat(list.w_rd[w].letter[l].segment[1]);
      list0.w_rd[w].letter[l].segment[1].velar:=velar(list.w_rd[w].letter[l].segment[1]);

      if list.w_rd[w].letter[l].segment[2].head='   '
        then begin
               list0.w_rd[w].letter[l].segment[2].sound:=list0.w_rd[w].letter[l].segment[1].sound;
               list0.w_rd[w].letter[l].segment[2].voice:=list0.w_rd[w].letter[l].segment[1].voice;
               list0.w_rd[w].letter[l].segment[2].apic :=list0.w_rd[w].letter[l].segment[1].apic;
               list0.w_rd[w].letter[l].segment[2].nasal:=list0.w_rd[w].letter[l].segment[1].nasal;
               list0.w_rd[w].letter[l].segment[2].palat:=list0.w_rd[w].letter[l].segment[1].palat;
               list0.w_rd[w].letter[l].segment[2].velar:=list0.w_rd[w].letter[l].segment[1].velar;
             end
        else begin
               list0.w_rd[w].letter[l].segment[2].sound:=sound(list.w_rd[w].letter[l].segment[2]);
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

procedure initcomp1;
var
  fp  : text;
  i,j : integer;
  ch  : char;
  v   : real;
begin
  openr(fp,cfile1);

  readln(fp,comp1.number);

  for i:=1 to comp1.number do begin
    initsegment(comp1.segment[i]);
    read(fp,ch);
    rsegment(fp,comp1.segment[i],ch);
    readln(fp);
    comp1.matrix[i,i]:=0;
  end;

  for i:=2 to comp1.number do begin
    for j:=1 to (i-1) do begin
      readln(fp,v);
      comp1.matrix[i,j]:=ln(v+1);
      comp1.matrix[j,i]:=ln(v+1);
    end;
  end;

  close(fp);
end;

procedure initcomp2;
var
  fp  : text;
  i,j : integer;
  ch  : char;
  v   : real;
begin
  openr(fp,cfile2);

  readln(fp,comp2.number);

  for i:=1 to comp2.number do begin
    initstress(comp2.stress[i]);

    if not eoln(fp)
      then begin
             read(fp,ch);
             rstress(fp,comp2.stress[i],ch);
           end
      else {nothing};

    readln(fp);
    comp2.matrix[i,i]:=0;
  end;

  for i:=2 to comp2.number do begin
    for j:=1 to (i-1) do begin
      readln(fp,v);
      comp2.matrix[i,j]:=v;
      comp2.matrix[j,i]:=v;
    end;
  end;

  close(fp);
end;

procedure initvc1;
var
  i,j : integer;
begin
  vc1:=0;
  for i:=2 to comp1.number do begin
    for j:=1 to (i-1) do begin
      if comp1.matrix[i,j]>vc1
        then vc1:=comp1.matrix[i,j]
        else {nothing};
    end;
  end;
end;

procedure initvc2;
var
  i,j : integer;
begin
  vc2:=0;
  for i:=2 to comp2.number do begin
    for j:=1 to (i-1) do begin
      if comp2.matrix[i,j]>vc2
        then vc2:=comp2.matrix[i,j]
        else {nothing};
    end;
  end;
end;

procedure initwl1;
var
  list     : listtype0;
  w,l,t,w0 : integer;
  fp       : text;
begin
  for w:=1 to m do begin
    for l:=0 to ml do begin
      wl1[w,l]:=0;
    end;
  end;
  
  for t:=1 to n do begin
    openr(fp,path+dialect[t]);

    lnr:=0;m:=0;
    while not eof(fp) do begin
      inc(lnr);
      rline(fp,list);

      if (f[lnr]<>0)
        then begin
	       m:=m+1;
	       
               if (list.number>0)
	         then for w0:=1 to list.number do begin
                        wl1[m,0]:=wl1[m,0]+(1/list.number);
                        wl1[m,list.w_rd[w0].lexeme]:=wl1[m,list.w_rd[w0].lexeme]+(1/list.number);
                      end
		 else {nothing};
	     end
	else {nothing};
    end;

    close(fp);
  end;

  for w:=1 to m do begin
    for l:=1 to ml do begin
      wl1[w,l]:=(wl1[w,l]/wl1[w,0])*100;
    end;
  end;
end;

procedure initwl2;
var
  list     : listtype0;
  w,l,t,w0 : integer;
  fp       : text;
begin
  for w:=1 to m do begin
    for l:=0 to ml do begin
      wl2[w,l]:=0;
    end;
  end;
  
  for t:=1 to n do begin
    openr(fp,path+dialect[t]);

    lnr:=0;m:=0;
    while not eof(fp) do begin
      inc(lnr);
      rline(fp,list);

      if (f[lnr]<>0)
        then begin
	       m:=m+1;
	       
               if (list.number>0)
	         then for w0:=1 to list.number do begin
                        wl2[m,0]:=wl2[m,0]+(1/list.number);
                        wl2[m,list.w_rd[w0].morpheme]:=wl2[m,list.w_rd[w0].morpheme]+(1/list.number);
                      end
		 else {nothing};
	     end
	else {nothing};
    end;

    close(fp);
  end;

  for w:=1 to m do begin
    for l:=1 to ml do begin
      wl2[w,l]:=(wl2[w,l]/wl2[w,0])*100;
    end;
  end;
end;

procedure determine1(word1,word2:wordtype0;var distance:real);
begin
  if (word1.lexeme=word2.lexeme)
    then case giw1 of
           1 : distance:=0;
           2 : distance:=(wl1[m,word1.lexeme]+wl1[m,word2.lexeme])/2;
	 end
    else distance:=100;
end;      

procedure determine2(word1,word2:wordtype0;var distance:real);
begin
  if (word1.lexeme=word2.lexeme)
    then if (word1.morpheme=word2.morpheme)
           then case giw2 of
                  1 : distance:=0;
                  2 : distance:=(wl2[m,word1.morpheme]+wl2[m,word2.morpheme])/2;
                end
           else distance:=100
    else distance:=101;
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

function sub10(segment1,segment2:segmenttype0):real;
var
  d : real;
  k : integer;
begin
  if (    vowel1(segment1.sound) and     vowel1(segment2.sound)) or
     (consonant1(segment1.sound) and consonant1(segment2.sound)) or
     (    vowel2(segment1.sound) and consonant2(segment2.sound)) or
     (consonant2(segment1.sound) and     vowel2(segment2.sound))
    then begin
           d:=comp1.matrix[segment1.sound,segment2.sound];
           k:=1;

           if (segment1.voice) and (not segment2.voice)
             then begin
                    d:=d+comp1.matrix[opposite(segment1.sound),segment2.sound];
                    k:=k+1;
                  end
             else
           if (not segment1.voice) and (segment2.voice)
             then begin
                    d:=d+comp1.matrix[segment1.sound,opposite(segment2.sound)];
                    k:=k+1;
                  end
             else
           if (segment1.voice) and (segment2.voice)
             then begin
                    d:=d+comp1.matrix[opposite(segment1.sound),opposite(segment2.sound)];
                    k:=k+1;
                  end
             else {nothing};

           if (segment1.apic) and (not segment2.apic)
             then if (segment1.sound=61)
                    then begin
                           d:=d+comp1.matrix[65,segment2.sound];
                           k:=k+1;
                         end
                    else
                  if (segment1.sound=62)
                    then begin
                           d:=d+comp1.matrix[66,segment2.sound];
                           k:=k+1;
                         end
                    else halt0('apical applied to wrong sound')
             else
           if (not segment1.apic) and (segment2.apic)
             then if (segment2.sound=61)
                    then begin
                           d:=d+comp1.matrix[segment1.sound,65];
                           k:=k+1;
                         end
                    else
                  if (segment2.sound=62)
                    then begin
                           d:=d+comp1.matrix[segment1.sound,66];
                           k:=k+1;
                         end
                    else halt0('apical applied to wrong sound')
             else
           if (segment1.apic) and (segment2.apic)
             then if (segment1.sound=61) and (segment2.sound=62)
                    then begin
                           d:=d+comp1.matrix[65,66];
                           k:=k+1;
                         end
                    else
                  if (segment1.sound=62) and (segment2.sound=61)
                    then begin
                           d:=d+comp1.matrix[66,65];
                           k:=k+1;
                         end
                    else
                  if (segment1.sound=61) and (segment2.sound=61)
                    then begin
                           d:=d+comp1.matrix[65,65];
                           k:=k+1;
                         end
                    else
                  if (segment1.sound=62) and (segment2.sound=62)
                    then begin
                           d:=d+comp1.matrix[66,66];
                           k:=k+1;
                         end
                    else halt0('apical applied to wrong sound')
             else {nothing};

           if (segment1.nasal) and (not segment2.nasal)
             then begin
                    d:=d+comp1.matrix[45,segment2.sound];
                    k:=k+1;
                  end
             else
           if (not segment1.nasal) and (segment2.nasal)
             then begin
                    d:=d+comp1.matrix[segment1.sound,45];
                    k:=k+1;
                  end
             else
           if (segment1.nasal) and (segment2.nasal)
             then begin
                    d:=d+(comp1.matrix[45,45]);
                    k:=k+1;
                  end
             else {nothing};

           if (segment1.palat) and (not segment2.palat)
             then begin
                    d:=d+comp1.matrix[82,segment2.sound];
                    k:=k+1;
                  end
             else
           if (not segment1.palat) and (segment2.palat)
             then begin
                    d:=d+comp1.matrix[segment1.sound,82];
                    k:=k+1;
                  end
             else
           if (segment1.palat) and (segment2.palat)
             then begin
                    d:=d+(comp1.matrix[82,82]);
                    k:=k+1;
                  end
             else {nothing};

           if (segment1.velar) and (not segment2.velar)
             then begin
                    d:=d+comp1.matrix[70,segment2.sound];
                    k:=k+1;
                  end
             else
           if (not segment1.velar) and (segment2.velar)
             then begin
                    d:=d+comp1.matrix[segment1.sound,70];
                    k:=k+1;
                  end
             else
           if (segment1.velar) and (segment2.velar)
             then begin
                    d:=d+(comp1.matrix[70,70]);
                    k:=k+1;
                  end
             else {nothing};

           sub10:=d/k;
         end
    else sub10:=maxint;
end;

function sub1(letter1,letter2:lettertype0):real;
begin
  sub1:=(sub10(letter1.segment[1],letter2.segment[1])+
         sub10(letter1.segment[2],letter2.segment[2]))/2;
end;

function sub2(letter1,letter2:lettertype0):real;
begin
  if vowel1(letter1.segment[1].sound) and vowel1(letter2.segment[1].sound) and
     vowel1(letter1.segment[2].sound) and vowel1(letter2.segment[2].sound)  
    then sub2:=comp2.matrix[letter1.stress,letter2.stress]
    else sub2:=0;
end;

function insdel10(segment:segmenttype0):real;
var
  d : real;
  k : integer;
begin
  d:=comp1.matrix[1,segment.sound];
  k:=1;

  if (segment.voice)
    then begin
           d:=d+comp1.matrix[1,opposite(segment.sound)];
           k:=k+1;
         end
    else {nothing};

 if (segment.apic)
   then if segment.sound=61
          then begin
                 d:=d+comp1.matrix[1,65];
                 k:=k+1;
               end
          else
        if segment.sound=62
          then begin
                 d:=d+comp1.matrix[1,66];
                 k:=k+1;
               end
          else halt0('apical applied to wrong sound')
   else {nothing};

  if (segment.nasal)
    then begin
           d:=d+comp1.matrix[1,45];
           k:=k+1;
         end
    else {nothing};

  if (segment.palat)
    then begin
           d:=d+comp1.matrix[1,82];
           k:=k+1;
         end
    else {nothing};

  if (segment.velar)
    then begin
           d:=d+comp1.matrix[1,70];
           k:=k+1;
         end
    else {nothing};

  insdel10:=d/k;
end;

function insdel1(letter:lettertype0):real;
begin
  insdel1:=(insdel10(letter.segment[1])+
            insdel10(letter.segment[2]))/2;
end;

function insdel2(letter:lettertype0):real;
begin
  insdel2:=0;
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

procedure check(var value1:real;upper1,upperleft1,left1:real;
                var value2:real;upper2,upperleft2,left2:real;
                var value3:real;upper3,upperleft3,left3:real);
begin
  value1:=min(upper1,upperleft1,left1);

  if upper1<>value1
    then upper2:= maxint;

  if upperleft1<>value1
    then upperleft2:= maxint;

  if left1<>value1
    then left2:= maxint;

  value2:=min(upper2,upperleft2,left2);

  if upper2<>value2
    then upper3:=-maxint;

  if upperleft2<>value2
    then upperleft3:=-maxint;

  if left2<>value2
    then left3:=-maxint;

  value3:=max(upper3,upperleft3,left3);
end;

procedure addtrace(upper1,upperleft1,left1,
                   upper2,upperleft2,left2,
                   upper3,upperleft3,left3:real;l1,l2:integer);
var
  pointer : integer;
begin
  pointer:=0;

  if (upper1=matrix1[l1,l2]) and (upper2=matrix2[l1,l2]) and (upper3=matrix3[l1,l2])
    then pointer:=pointer+8;

  if (upperleft1=matrix1[l1,l2]) and (upperleft2=matrix2[l1,l2]) and (upperleft3=matrix3[l1,l2])
    then pointer:=pointer+4;

  if (left1=matrix1[l1,l2]) and (left2=matrix2[l1,l2]) and (left3=matrix3[l1,l2])
    then pointer:=pointer+2;

  trace[l1,l2]:=pointer;
end;

procedure addalign(l1,l2:integer;d:real;a:integer);
begin
  align[a].word1:=l1;
  align[a].word2:=l2;
  align[a].dist :=d ;
end;

procedure procalign(a:integer;word1,word2:wordtype0;var distance1,distance3:real);
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

  distance1:=(distance10*100)/vc1;
  distance3:=l0;  
end;

procedure proctrace(l1,l2,a:integer;word1,word2:wordtype0;var distance1,distance3:real);
var
  pointer : integer;
begin
  if (l1>0) or (l2>0)
    then begin
           pointer:=trace[l1,l2];

           if ((pointer= 8) or (pointer=10) or (pointer=12) or (pointer=14)) and (l1>0)
             then begin
                    addalign(l1,0,matrix1[l1,l2]-matrix1[l1-1,l2],a);
                    proctrace(l1-1,l2,a+1,word1,word2,distance1,distance3);
                    pointer:=pointer-8;
                  end;

           if ((pointer= 4) or (pointer= 6) or (pointer=12) or (pointer=14)) and (l1>0) and (l2>0)
             then begin
                    addalign(l1,l2,matrix1[l1,l2]-matrix1[l1-1,l2-1],a);
                    proctrace(l1-1,l2-1,a+1,word1,word2,distance1,distance3);
                    pointer:=pointer-4;
                  end;

           if ((pointer= 2) or (pointer= 6) or (pointer=10) or (pointer=14)) and (l2>0)
             then begin
                    addalign(0,l2,matrix1[l1,l2]-matrix1[l1,l2-1],a);
                    proctrace(l1,l2-1,a+1,word1,word2,distance1,distance3);
                    pointer:=pointer-2;
                  end;
         end
    else procalign(a,word1,word2,distance1,distance3);
end;

procedure determine3(word1,word2:wordtype0;var distance1,distance2,distance3:real);
var
  l1,l2                   : integer;
  left1,upperleft1,upper1 : real;
  left2,upperleft2,upper2 : real;
  left3,upperleft3,upper3 : real;
  m1,m2,m3                : real;
begin
  if (word1.lexeme=word2.lexeme) and (word1.morpheme=word2.morpheme)
    then begin
           for l1:=0 to word1.number do begin
             for l2:=0 to word2.number do begin
               upper1:= maxint;upperleft1:= maxint;left1:= maxint;
               upper2:= maxint;upperleft2:= maxint;left2:= maxint;
               upper3:=-maxint;upperleft3:=-maxint;left3:=-maxint;

               {deletion of word1[l1]}
               if (l1>0)
                 then begin
                        upper1:=matrix1[l1-1,l2]+insdel1(word1.letter[l1]);
                        upper2:=matrix2[l1-1,l2]+insdel2(word1.letter[l1]);
                        upper3:=matrix3[l1-1,l2]+1;
                      end;

               {substitution of word1[l1] by word2[l2]}
               if (l1>0) and (l2>0)
                 then begin
                        upperleft1:=matrix1[l1-1,l2-1]+sub1(word1.letter[l1],word2.letter[l2]);
                        upperleft2:=matrix2[l1-1,l2-1]+sub2(word1.letter[l1],word2.letter[l2]);
                        upperleft3:=matrix3[l1-1,l2-1]+1;
                      end;

               {insertion of word2[l2]}
               if (l2>0)
                 then begin
                        left1:=matrix1[l1,l2-1]+insdel1(word2.letter[l2]);
                        left2:=matrix2[l1,l2-1]+insdel2(word2.letter[l2]);
                        left3:=matrix3[l1,l2-1]+1;
                      end;

               check(matrix1[l1,l2],upper1,upperleft1,left1,
                     matrix2[l1,l2],upper2,upperleft2,left2,
                     matrix3[l1,l2],upper3,upperleft3,left3);

               if matrix1[l1,l2]= maxint then matrix1[l1,l2]:=0;
               if matrix2[l1,l2]= maxint then matrix2[l1,l2]:=0;
               if matrix3[l1,l2]=-maxint then matrix3[l1,l2]:=0;

               addtrace(upper1,upperleft1,left1,
                        upper2,upperleft2,left2,
                        upper3,upperleft3,left3,l1,l2);
             end;
           end;

           m1:=matrix1[word1.number,word2.number];
           m2:=matrix2[word1.number,word2.number];
           m3:=matrix3[word1.number,word2.number];

           distance1:=(m1*100)/vc1;
           distance2:=(m2*100)/vc2;
           distance3:=m3;	   

           proctrace(word1.number,word2.number,1,word1,word2,distance1,distance3);
         end
    else distance1:=20001;
end;

procedure calculate(i,j:integer;var d_lex,d_mor,d_pho,d_ton,d_len:dtype);
var
  list1,list2                   : listtype0;
  w1,w2                         : integer;
  dialect1,dialect2             : text;
  matrix1,matrix2,matrix3       : array[1..10] of array[1..10] of real;
  c1,c2                         : array[0..10] of integer;
  h1,h2                         : integer;
  d10,d20,d30                   : real;
  w10,w20                       : integer;
  n_lex,n_mor,n_pho,n_ton,n_len : integer;
begin
  openr(dialect1,path+dialect[i]);
  openr(dialect2,path+dialect[j]);

  lnr:=0; m:=0; m0_lex:=0; m0_mor:=0; m0_pho:=0; m0_ton:=0; m0_len:=0; m00:=0;
  while (not eof(dialect1)) and (not eof(dialect2)) do begin
    inc(lnr);

    rline(dialect1,list1);
    rline(dialect2,list2);

    if (f[lnr]<>0)
      then begin
             m:=m+1;
             if (list1.number>0) and (list2.number>0)
               then begin
                      {calculate lexical distance}

                      for w1:=1 to list1.number do
                        for w2:=1 to list2.number do
                          determine1(list1.w_rd[w1],list2.w_rd[w2],matrix1[w1,w2]);

                      for w1:=1 to list1.number do
                        c1[w1]:=list2.number;
                      for w2:=1 to list2.number do
                        c2[w2]:=list1.number;

                      d_lex[m]:=0;
                      n_lex:=0;
                      
                      for h1:=1 to list1.number do begin
                        for h2:=1 to list2.number do begin

                          d10:=1.7e+38;
                          for w1:=1 to list1.number do begin
                            for w2:=1 to list2.number do begin
                              if (matrix1[w1,w2]<d10) and (c1[w1]>0) and (c2[w2]>0)
                                then begin
                                       d10:=matrix1[w1,w2];
                                       w10:=w1;
                                       w20:=w2;
                                     end
                                else {nothing};
                            end;
                          end;

                          d_lex[m]:=d_lex[m]+d10;
                          n_lex:=n_lex+1;

                          c1[w10]:=c1[w10]-1;
                          c2[w20]:=c2[w20]-1;
                        end;
                      end;

                      d_lex[m]:=(d_lex[m]/n_lex)*f[lnr]; 
                      m0_lex:=m0_lex+f[lnr];

                      {calculate morphological distance}

                      for w1:=1 to list1.number do
                        for w2:=1 to list2.number do
                          determine2(list1.w_rd[w1],list2.w_rd[w2],matrix1[w1,w2]);

                      for w1:=1 to list1.number do
                        c1[w1]:=list2.number;
                      for w2:=1 to list2.number do
                        c2[w2]:=list1.number;

                      d_mor[m]:=0;
                      n_mor:=0;
                      
                      for h1:=1 to list1.number do begin
                        for h2:=1 to list2.number do begin

                          d10:=1.7e+38;
                          for w1:=1 to list1.number do begin
                            for w2:=1 to list2.number do begin
                              if (matrix1[w1,w2]<d10) and (c1[w1]>0) and (c2[w2]>0)
                                then begin
                                       d10:=matrix1[w1,w2];
                                       w10:=w1;
                                       w20:=w2;
                                     end
                                else {nothing};
                            end;
                          end;

                          if (d10<101)
                            then begin
                                   d_mor[m]:=d_mor[m]+d10;
                                   n_mor:=n_mor+1;
                                 end
                            else {nothing};

                          c1[w10]:=c1[w10]-1;
                          c2[w20]:=c2[w20]-1;
                        end;
                      end;

                      if (n_mor>0)
                        then begin
                               d_mor[m]:=(d_mor[m]/n_mor)*f[lnr];
                               m0_mor:=m0_mor+f[lnr];
                             end
                        else begin
                               d_mor[m]:=-1*f[lnr];
                             end;

                      {calculate pronunciation distance}

                      for w1:=1 to list1.number do
                        for w2:=1 to list2.number do
                          determine3(list1.w_rd[w1],list2.w_rd[w2],matrix1[w1,w2],matrix2[w1,w2],matrix3[w1,w2]);

                      for w1:=1 to list1.number do
                        c1[w1]:=list2.number;
                      for w2:=1 to list2.number do
                        c2[w2]:=list1.number;

                      d_pho[m]:=0;n_pho:=0;
                      d_ton[m]:=0;n_ton:=0;
                      d_len[m]:=0;n_len:=0;

                      for h1:=1 to list1.number do begin
                        for h2:=1 to list2.number do begin

                          d10:=1.7e+38;
                          for w1:=1 to list1.number do begin
                            for w2:=1 to list2.number do begin
                              if (matrix1[w1,w2]<d10) and (c1[w1]>0) and (c2[w2]>0)
                                then begin
                                       d10:=matrix1[w1,w2];
                                       d20:=matrix2[w1,w2];
                                       d30:=matrix3[w1,w2];
                                       w10:=w1;
                                       w20:=w2;
                                     end
                                else

                              if (matrix1[w1,w2]=d10) and (c1[w1]>0) and (c2[w2]>0)
                                 and (matrix2[w1,w2]<matrix2[w10,w20])
                                then begin
                                       d20:=matrix2[w1,w2];
                                       d30:=matrix3[w1,w2];
                                       w10:=w1;
                                       w20:=w2;
                                     end
                                else {nothing};
                            end;
                          end;

                          if (d10<20001)
                            then begin
                                   d_pho[m]:=d_pho[m]+d10;
                                   d_ton[m]:=d_ton[m]+d20;
                                   d_len[m]:=d_len[m]+d30;

                                   n_pho:=n_pho+1;
                                   n_ton:=n_ton+1;
                                   n_len:=n_len+1;
                                 end
                            else {nothing};

                          c1[w10]:=c1[w10]-1;
                          c2[w20]:=c2[w20]-1;
                        end;
                      end;

                      if (n_pho>0)
                        then begin
                               d_pho[m]:=(d_pho[m]/n_pho)*f[lnr];
                               d_ton[m]:=(d_ton[m]/n_ton)*f[lnr];
                               d_len[m]:=(d_len[m]/n_len)*f[lnr];

                               m0_pho:=m0_pho+f[lnr];
                               m0_ton:=m0_ton+f[lnr];
                               m0_len:=m0_len+f[lnr];
                             end
                        else begin
                               d_pho[m]:=-1*f[lnr];
                               d_ton[m]:=-1*f[lnr];
                               d_len[m]:=-1*f[lnr];
                             end;
                    end
               else begin
                      d_lex[m]:=-1*f[lnr];
                      d_mor[m]:=-1*f[lnr];
                      d_pho[m]:=-1*f[lnr];
                      d_ton[m]:=-1*f[lnr];
                      d_len[m]:=-1*f[lnr];      
                    end;
             m00:=m00+f[lnr];
           end
      else {nothing}
  end;

  close(dialect1);
  close(dialect2);
end;

procedure process1(var fp:text;d:dtype;m0:real);
var
  w        : integer;
  sum,mean : real;
begin
  sum:=0;

  for w:=1 to m do begin
    if d[w]>=0
      then sum:=sum+d[w]
      else {nothing}
  end;

  mean:=sum/m0;

  for w:=1 to m do begin
    if d[w]<0
      then d[w]:=d[w]*mean
      else {nothing};
  end;

  for w:=1 to m do begin
    writeln(fp,d[w]/m00);
  end;
end;

procedure process2(var fp:text;d:dtype;m0:real;d_len:dtype;m0_len:real);
var
  w                         : integer;
  sum,sum_len,mean,mean_len : real;
begin
  sum:=0;

  for w:=1 to m do begin
    if d[w]>=0
      then sum:=sum+d[w]
      else {nothing}
  end;

  mean:=sum/m0;

  for w:=1 to m do begin
    if d[w]<0
      then d[w]:=d[w]*mean
      else {nothing};
  end;

  sum_len:=0;

  for w:=1 to m do begin
    if d_len[w]>=0
      then sum_len:=sum_len+d_len[w]
      else {nothing};
  end;

  mean_len:=sum_len/m0_len;

  for w:=1 to m do begin
    if d_len[w]<0
      then d_len[w]:=-1*d_len[w]*mean_len
      else {nothing};
  end;

  sum_len:=0;

  for w:=1 to m do begin
    if d_len[w]>=0
      then sum_len:=sum_len+d_len[w]
      else halt;
  end;

  for w:=1 to m do begin
    writeln(fp,d[w]/sum_len);
  end;
end;

begin {main}
  getparameters;

  initpath;
  initdialect;
  initfreq;

  initcomp1;
  initcomp2;

  initvc1;
  initvc2;

  initwl1;
  initwl2;

  openw(fp_lex,'data.lex.dis');
  openw(fp_mor,'data.mor.dis');
  openw(fp_pho,'data.pho.dis');
  openw(fp_ton,'data.ton.dis');

  writeln(fp_lex,n);
  writeln(fp_mor,n);
  writeln(fp_pho,n);
  writeln(fp_ton,n);

  writeln(fp_lex,m);
  writeln(fp_mor,m);
  writeln(fp_pho,m);
  writeln(fp_ton,m);

  writeln(stderr);
  for i:=2 to n do begin
    for j:=1 to (i-1) do begin
      writeln(stderr,'Processing ',dialect[i],' and ',dialect[j]);
      calculate(i,j,d_lex,d_mor,d_pho,d_ton,d_len);

      process1(fp_lex,d_lex,m0_lex);
      process1(fp_mor,d_mor,m0_mor);
      process2(fp_pho,d_pho,m0_pho,d_len,m0_len);
      process2(fp_ton,d_ton,m0_ton,d_len,m0_len);

      writeln(output,m0_lex,' ',m0_mor,'  ',m0_pho,' ',m0_ton);
    end;
    writeln(stderr);
  end;

  close(fp_lex);
  close(fp_mor);
  close(fp_pho);
  close(fp_ton);
end   {main}.
