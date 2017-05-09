Unit MainForm;

Interface

Uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants,
  System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.ExtCtrls, Vcl.StdCtrls,
  Vcl.ComCtrls, Generics.Collections,
  VCFRecord;

Type
  EVCFIncidentType = (
    vcfitTruePositive,
    vcfitFalseNegative,
    vcfitFalsePositive
  );
  EVCFIncidentTypeSet = Set Of EVCFIncidentType;

  TVCFIncident = Record
    IncidentType : EVCFIncidentType;
    VCFRecord : TVCFRecord;
    end;

  TMainFrm = class(TForm)
    MainPanel: TPanel;
    SetsAndFilesGroupBox: TGroupBox;
    SettingsGroupBox: TGroupBox;
    ResultListView: TListView;
    TruthSetFileEdit: TEdit;
    TestSetFileEdit: TEdit;
    Label1: TLabel;
    Label2: TLabel;
    TruthSetBrowseButton: TButton;
    TestSetBrowseButton: TButton;
    SNPsCheckBox: TCheckBox;
    InsertionsCheckBox: TCheckBox;
    DeletionsCheckBox: TCheckBox;
    TruePositivesCheckBox: TCheckBox;
    FalseNegativesCheckBox: TCheckBox;
    FalsePositivesCheckBox: TCheckBox;
    TruthSetOpenDialog: TOpenDialog;
    TestSetOpenDialog: TOpenDialog;
    procedure FormCreate(Sender: TObject);
    procedure FormDestroy(Sender: TObject);
    procedure FilterCheckBoxClick(Sender: TObject);
  Private
    FTruthFile : TVCFFile;
    FTestFile : TVCFFile;
    FFilteredTruth : TList<TVCFRecord>;
    FFilteredTest : TList<TVCFRecord>;
    FFilterSet : EVCFRecordTypeSet;
    FFilterSetIncidents : EVCFIncidentTypeSet;
    FIncidents : TList<TVCFIncident>;
    Procedure BuildFilterSet;
    Procedure ApplyFilterSet(AFile:TVCFFile; AList:TList<TVCFRecord>);
    Procedure ProcessIncidents;
    Procedure AddIncident(AType:EVCFIncidentType; ARecord:TVCFRecord; AList:TList<TVCFIncident>);
  end;

Var
  MainFrm: TMainFrm;

Implementation

Procedure TMainFrm.BuildFilterSet;
begin
FFilterSet := [];
If SNPsCheckBox.Checked Then
  FFilterSet := FFilterSet + [vcfrtSNP];

If DeletionsCheckBox.Checked Then
  FFilterSet := FFilterSet + [vcfrtDeletion];

If InsertionsCheckBox.Checked Then
  FFilterSet := FFilterSet + [vcfrtInsertion];

FFilterSetIncidents := [];
If TruePositivesCheckBox.Checked Then
  FFilterSetIncidents := FFilterSetIncidents + [vcfitTruePositive];

If FalseNegativesCheckBox.Checked Then
  FFilterSetIncidents := FFilterSetIncidents + [vcfitFalseNegative];

If FalsePositivesCheckBox.Checked Then
  FFilterSetIncidents := FFilterSetIncidents + [vcfitFalsePositive];
end;

Procedure TMainFrm.ApplyFilterSet(AFile:TVCFFile; AList:TList<TVCFRecord>);
begin
If Assigned(AFile) Then
  begin
  AList.Clear;
  AFile.Filter(FFilterSet, AList);
  end;
end;

Procedure TMainFrm.AddIncident(AType:EVCFIncidentType; ARecord:TVCFRecord; AList:TList<TVCFIncident>);
Var
  incident : TVCFIncident;
begin
If (AType In FFilterSetIncidents) Then
  begin
  incident.IncidentType := AType;
  incident.VCFRecord := ARecord;
  AList.Add(incident);
  end;
end;

Procedure TMainFrm.ProcessIncidents;
Var
  incident : TVCFIncident;
  r1, r2 : TVCFRecord;
  index1, index2 : Integer;
  newList : TList<TVCFIncident>;
begin
newList := TList<TVCFIncident>.Create;
index1 := 0;
index2 := 0;
While (index1 < FFilteredTruth.Count) And (index2 < FFilteredTest.Count) Do
  begin
  r1 := FFilteredTruth[index1];
  r2 := FFilteredTest[index2];
  If TVCFRecord.Equal(r1, r2) Then
    begin
    AddIncident(vcfitTruePositive, r2, newList);
    Inc(index1);
    Inc(index2);
    end
  Else If r1.Pos < r2.Pos Then
    begin
    AddIncident(vcfitFalseNegative, r1, newList);
    Inc(index1);
    end
  Else If r1.Pos > r2.Pos Then
    begin
    AddIncident(vcfitFalsePositive, r2, newList);
    Inc(index2);
    end
  Else begin
    AddIncident(vcfitFalseNegative, r1, newList);
    AddIncident(vcfitFalsePositive, r2, newList);
    Inc(index1);
    Inc(index2);
    end;
  end;

While index1 < FFilteredTruth.Count Do
  begin
  r1 := FFilteredTruth[index1];
  AddIncident(vcfitFalseNegative, r1, newList);
  Inc(index1);
  end;

While index2 < FFilteredTest.Count Do
  begin
  r2 := FFilteredTest[index2];
  AddIncident(vcfitFalsePositive, r2, newList);
  Inc(index2);
  end;

newList.Free;
end;

{$R *.DFM}

Procedure TMainFrm.FilterCheckBoxClick(Sender: TObject);
begin
BuildFilterSet;
ApplyFilterSet(FTruthFile, FFilteredTruth);
ApplyFilterSet(FTestFile, FFilteredTest);
ProcessIncidents;
end;

Procedure TMainFrm.FormCreate(Sender: TObject);
begin
FFIlterSet := [];
FFilteredTruth := TList<TVCFRecord>.Create;
FFilteredTest := TList<TVCFRecord>.Create;
FIncidents := TList<TVCFIncident>.Create;
BuildFilterSet;
end;

Procedure TMainFrm.FormDestroy(Sender: TObject);
begin
ResultListView.Items.Count := 0;
FIncidents.Free;
FFilteredTest.Free;
FFilteredTruth.Free;
end;

End.

