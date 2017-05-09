Unit VCFRecord;

Interface

Uses
  Classes, Generics.Collections;

Type
  EVCFRecordType = (
    vcfrtSNP,
    vcfrtInsertion,
    vcfrtDeletion,
    vcfrtReplace
  );

   EVCFRecordTypeSet = Set Of EVCFRecordType;

  TVCFRecord = Class
    Private
      FChrom : WideString;
      FPos : UInt64;
      FId : WideString;
      FRef : WideString;
      FAlt : WideString;
      FInfo : WideString;
      FFilter : WideString;
      FFormat : WideString;
      FSample : WideString;

      FRW : Cardinal;
      FRC : Cardinal;
      FAW : Cardinal;
      FAC : Cardinal;

      FType : EVCFRecordType;
      Procedure SplitLine(ALine:WideString; ADelimiter:WideChar; AList:TStrings);
      Procedure ParseInfo;
    Public
      Class Function Equal(A:TVCFRecord; B:TVCFRecord):Boolean;

      Constructor Create(ALine:WideString); Reintroduce;


      Property Chrom : WideString Read FChrom;
      Property Pos : UInt64 Read FPos;
      Property Id : WideString Read FId;
      Property Ref : WideString Read FRef;
      Property Alt : WideString Read FAlt;
      Property Info : WideString Read FInfo;
      Property Filter : WideString Read FFilter;
      Property Format : WideString Read FFormat;
      Property Sample : WideString Read FSample;

      Property RecordType : EVCFRecordType Read FType;
      Property RW : Cardinal Read FRW;
      Property RC : Cardinal Read FRC;
      Property AW : Cardinal Read FAW;
      Property AC : Cardinal Read FAC;
    end;

    TVCFFile = Class
      Private
        FRecords : TObjectList<TVCFRecord>;
      Protected
        Function GetRecord(AIndex:Integer):TVCFRecord;
        Function GetRecordCount:Integer;
      Public
        Constructor Creae(AFileName:WideString); Reintroduce;
        Destructor Destroy; Override;

        Procedure Filter(AFilterSet:EVCFRecordTypeSet; AResult:TList<TVCFRecord>);

        Property Records [Index:integer] : TVCFRecord Read GetRecord;
        Property Count : Integer Read GetRecordCount;
      end;


Implementation

Uses
  SysUtils;

(** TVCFRecord **)

Procedure TVCFRecord.ParseInfo;
Var
  name : WideString;
  value : WideString;
  index : Integer;
  t : WideString;
  tuples : TStringList;
begin
tuples := TStringList.Create;
SplitLine(FInfo, ';', tuples);
For t In tuples DO
  begin
  index := System.Pos('=', t);
  name := Copy(t, 1, index - 1);
  value := Copy(t, index + 1, Length(t) - index);
  If name = 'RW' Then
    FRW := StrToInt64(value)
  Else If name = 'RC' Then
    FRC := StrToInt64(value)
  Else If name = 'AW' Then
    FAW := StrToInt64(value)
  Else If name = 'AC' Then
    FAC := Trunc(StrToFloat(value));

  end;

tuples.Free;
end;

Procedure TVCFRecord.SplitLine(ALine:WideString; ADelimiter:WideChar; AList:TStrings);
Var
  tabIndex : Integer;
begin
Repeat
tabIndex := System.Pos(ADelimiter, ALine);
If tabIndex > 0 Then
  begin
  AList.Add(Copy(ALine, 1, tabIndex - 1));
  System.Delete(ALine, 1, tabIndex);
  end;
Until tabIndex <= 0;

If ALine <> '' Then
  AList.Add(ALine);
end;

Class Function TVCFRecord.Equal(A:TVCFRecord; B:TVCFRecord):Boolean;
begin
Result := (
  (A.FPos = B.FPos) And
  (A.FRef = B.FRef) And
  (A.FAlt = B.FAlt)
);
end;


Constructor TVCFRecord.Create(ALine:WideString);
Var
  refLen : Integer;
  altLen : Integer;
  fieldList : TStringList;
begin
Inherited Create;
If Length(ALine) = 0 Then
  Raise Exception.Create('Empty VCF record line');

If ALine[1] = '#' Then
  Raise Exception.Create('Header VCF line');

fieldList := TStringList.Create;
Try
  SplitLine(ALine, #9, fieldList);
  FChrom := fieldList[0];
  FPos := StrToInt64(fieldList[1]);
  FId := fieldList[2];
  FRef := fieldList[3];
  FAlt := fieldList[4];

  FFilter := fieldList[6];
  FInfo := fieldList[7];
  FFormat := fieldList[8];
  FSample := fieldList[9];
  refLen := Length(FRef);
  altLen := Length(FAlt);
  If (refLen = 1) And (altLen = 1) Then
    FType := vcfrtSNP
  Else If (refLen = 1) Then
    FType := vcfrtInsertion
  Else If (altLen = 1) Then
    FType := vcfrtDeletion
  Else FType := vcfrtReplace;

  ParseInfo;
Finally
  fieldList.Free;
  end;
end;


(** TVCFFile **)

Function TVCFFile.GetRecord(AIndex:Integer):TVCFRecord;
begin
Result := FRecords[AIndex];
end;

Constructor TVCFFile.Creae(AFileName:WideString);
Var
  I : Integer;
  notHeader : Boolean;
  lineIndex : Integer;
  line : WideString;
  lines : TStringList;
begin
Inherited Create;
FRecords := TObjectList<TVCFRecord>.Create;
lines := TStringList.Create;
lines.LoadFromFile(AFileName);
lineIndex := 0;
If lines.Count > 0 Then
  begin
  notHeader := False;
  Repeat
  line := lines[lineIndex];
  Inc(lineIndex);
  notHeader := (Length(line) > 0) And (line[1] <> '#');
  Until notHeader;

  For I := lineIndex - 1 To lines.Count - 1 Do
    begin
    line := lines[I];
    If Length(line) = 0 Then
      Continue;

    FRecords.Add(TVCFRecord.Create(line));
    end;
  end;

lines.Free;
end;

Destructor TVCFFile.Destroy;
begin
FRecords.Free;
Inherited Destroy;
end;

Function TVCFFile.GetRecordCount:Integer;
begin
Result := FRecords.Count;
end;


Procedure TVCFFile.Filter(AFilterSet:EVCFRecordTypeSet; AResult:TList<TVCFRecord>);
Var
  r : TVCFRecord;
begin
For r In FRecords Do
  begin
  If (r.RecordType In AFilterSet) Then
    AResult.Add(r);
  end;
end;

End.
