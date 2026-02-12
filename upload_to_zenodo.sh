#!/bin/bash

# ============================================================
# Zenodo Upload Script for Large Files
# ============================================================

# 配置（替换成你的 token）
ACCESS_TOKEN="3i9eh4KUOxi9I6EPbVlAlj216UKOWs5BnYi5qhNEmhxqx7Z62rfEgsHh0WDM"
FILE_PATH="data/example_data/combined_obj.RData"
FILE_NAME="combined_obj.RData"

# 元数据
TITLE="TLS Spatial Transcriptomics Full Dataset"
DESCRIPTION="Complete Visium spatial transcriptomics data for TLS analysis across multiple cancer types"
CREATORS='[{"name": "Ange.Y", "affiliation": "University of Tokyo"}]'
UPLOAD_TYPE="dataset"

echo "========================================="
echo "Zenodo Upload Script (Fixed)"
echo "========================================="
echo "File: $FILE_PATH"
echo "Size: $(du -h $FILE_PATH | cut -f1)"
echo ""

# 检查 jq 是否安装
if ! command -v jq &> /dev/null; then
    echo "❌ Error: jq is not installed"
    echo "Please install: sudo apt-get install jq"
    exit 1
fi

# 1. 创建新的 deposition
echo "[1/5] Creating new deposition..."
RESPONSE=$(curl -s -X POST \
  "https://zenodo.org/api/deposit/depositions" \
  -H "Content-Type: application/json" \
  -d '{}' \
  -H "Authorization: Bearer $ACCESS_TOKEN")

# 使用 jq 解析 JSON
DEPOSITION_ID=$(echo $RESPONSE | jq -r '.id')
BUCKET_URL=$(echo $RESPONSE | jq -r '.links.bucket')

if [ "$DEPOSITION_ID" == "null" ] || [ -z "$DEPOSITION_ID" ]; then
    echo "❌ Error: Failed to create deposition"
    echo "Response: $RESPONSE"
    exit 1
fi

echo "✅ Deposition created: ID = $DEPOSITION_ID"
echo "✅ Bucket URL: $BUCKET_URL"
echo ""

# 2. 上传文件
echo "[2/5] Uploading file..."
echo "⏳ This will take 10-30 minutes for 4.3GB file"
echo "Progress:"

# 使用 curl 的进度条
curl -X PUT \
  "$BUCKET_URL/$FILE_NAME" \
  -H "Authorization: Bearer $ACCESS_TOKEN" \
  -H "Content-Type: application/octet-stream" \
  --data-binary @"$FILE_PATH" \
  --progress-bar -o /dev/null

UPLOAD_STATUS=$?

if [ $UPLOAD_STATUS -ne 0 ]; then
    echo ""
    echo "❌ Upload failed with status: $UPLOAD_STATUS"
    exit 1
fi

echo ""
echo "✅ File uploaded successfully"
echo ""

# 3. 添加元数据
echo "[3/5] Adding metadata..."
METADATA=$(cat <<EOF
{
  "metadata": {
    "title": "$TITLE",
    "upload_type": "$UPLOAD_TYPE",
    "description": "$DESCRIPTION",
    "creators": $CREATORS,
    "access_right": "open",
    "license": "CC-BY-4.0",
    "keywords": ["spatial transcriptomics", "tertiary lymphoid structures", "TLS", "cancer immunology", "Visium", "EcoTyper"],
    "related_identifiers": [
      {
        "identifier": "https://github.com/anka2029/TLS_spatial_analysis",
        "relation": "isSupplementTo",
        "resource_type": "software"
      }
    ]
  }
}
EOF
)

METADATA_RESPONSE=$(curl -s -X PUT \
  "https://zenodo.org/api/deposit/depositions/$DEPOSITION_ID" \
  -H "Content-Type: application/json" \
  -H "Authorization: Bearer $ACCESS_TOKEN" \
  -d "$METADATA")

echo "✅ Metadata added"
echo ""

# 4. 获取 DOI
echo "[4/5] Getting pre-assigned DOI..."
DOI_RESPONSE=$(curl -s -X GET \
  "https://zenodo.org/api/deposit/depositions/$DEPOSITION_ID" \
  -H "Authorization: Bearer $ACCESS_TOKEN")

PRERESERVE_DOI=$(echo $DOI_RESPONSE | jq -r '.metadata.prereserve_doi.doi')
CONCEPT_DOI=$(echo $DOI_RESPONSE | jq -r '.conceptdoi // empty')

echo "✅ Pre-reserved DOI: $PRERESERVE_DOI"
if [ ! -z "$CONCEPT_DOI" ]; then
    echo "✅ Concept DOI: $CONCEPT_DOI"
fi
echo ""

# 5. 完成
echo "[5/5] Upload complete!"
echo ""
echo "========================================="
echo "📦 UPLOAD SUMMARY"
echo "========================================="
echo "Deposition ID: $DEPOSITION_ID"
echo "Pre-reserved DOI: $PRERESERVE_DOI"
echo "Preview URL: https://zenodo.org/deposit/$DEPOSITION_ID"
echo "File: $FILE_NAME ($(du -h $FILE_PATH | cut -f1))"
echo ""
echo "⚠️  IMPORTANT: The upload is NOT published yet!"
echo ""
echo "Next steps:"
echo "1. Visit: https://zenodo.org/deposit/$DEPOSITION_ID"
echo "2. Review the metadata and file"
echo "3. Click 'Publish' to make it public"
echo ""
echo "After publishing, update your GitHub README.md with:"
echo "[![DOI](https://zenodo.org/badge/DOI/$PRERESERVE_DOI.svg)](https://doi.org/$PRERESERVE_DOI)"
echo ""
echo "========================================="
echo "========================================="
echo "Zenodo Upload Script"
echo "========================================="
echo "File: $FILE_PATH"
echo "Size: $(du -h $FILE_PATH | cut -f1)"
echo ""

# 1. 创建新的 deposition（空的上传记录）
echo "[1/5] Creating new deposition..."
RESPONSE=$(curl -s -X POST \
  "https://zenodo.org/api/deposit/depositions" \
  -H "Content-Type: application/json" \
  -d '{}' \
  -H "Authorization: Bearer $ACCESS_TOKEN")

# 提取 deposition ID 和 bucket URL
DEPOSITION_ID=$(echo $RESPONSE | grep -o '"id":[0-9]*' | head -1 | grep -o '[0-9]*')
BUCKET_URL=$(echo $RESPONSE | grep -o '"bucket":"[^"]*' | grep -o 'http[^"]*')

if [ -z "$DEPOSITION_ID" ]; then
    echo "❌ Error: Failed to create deposition"
    echo "Response: $RESPONSE"
    exit 1
fi

echo "✅ Deposition created: ID = $DEPOSITION_ID"
echo ""

# 2. 上传文件（使用 bucket API，支持大文件）
echo "[2/5] Uploading file (this may take 10-30 minutes for 4.3GB)..."
echo "Progress will be shown below:"
curl -X PUT \
  "$BUCKET_URL/$FILE_NAME" \
  -H "Authorization: Bearer $ACCESS_TOKEN" \
  -H "Content-Type: application/octet-stream" \
  --data-binary @"$FILE_PATH" \
  --progress-bar | cat

echo ""
echo "✅ File uploaded successfully"
echo ""

# 3. 添加元数据
echo "[3/5] Adding metadata..."
METADATA=$(cat <<EOF
{
  "metadata": {
    "title": "$TITLE",
    "upload_type": "$UPLOAD_TYPE",
    "description": "$DESCRIPTION",
    "creators": $CREATORS,
    "access_right": "open",
    "license": "CC-BY-4.0",
    "keywords": ["spatial transcriptomics", "tertiary lymphoid structures", "cancer immunology", "Visium"]
  }
}
EOF
)

curl -s -X PUT \
  "https://zenodo.org/api/deposit/depositions/$DEPOSITION_ID" \
  -H "Content-Type: application/json" \
  -H "Authorization: Bearer $ACCESS_TOKEN" \
  -d "$METADATA" > /dev/null

echo "✅ Metadata added"
echo ""

# 4. 获取 DOI（预览）
echo "[4/5] Getting pre-assigned DOI..."
DOI_RESPONSE=$(curl -s -X GET \
  "https://zenodo.org/api/deposit/depositions/$DEPOSITION_ID" \
  -H "Authorization: Bearer $ACCESS_TOKEN")

DOI=$(echo $DOI_RESPONSE | grep -o '"doi":"[^"]*' | head -1 | cut -d'"' -f4)
PRERESERVE_DOI=$(echo $DOI_RESPONSE | grep -o '"conceptdoi":"[^"]*' | head -1 | cut -d'"' -f4)

echo "✅ Pre-reserved DOI: $PRERESERVE_DOI"
echo ""

# 5. 提供发布选项（不自动发布，让你检查）
echo "[5/5] Upload complete!"
echo ""
echo "========================================="
echo "📦 UPLOAD SUMMARY"
echo "========================================="
echo "Deposition ID: $DEPOSITION_ID"
echo "Pre-reserved DOI: $PRERESERVE_DOI"
echo "Preview URL: https://zenodo.org/deposit/$DEPOSITION_ID"
echo ""
echo "⚠️  IMPORTANT: The upload is NOT published yet!"
echo ""
echo "Next steps:"
echo "1. Visit: https://zenodo.org/deposit/$DEPOSITION_ID"
echo "2. Review the metadata and file"
echo "3. Click 'Publish' to make it public and activate the DOI"
echo ""
echo "To publish via command line (optional):"
echo "curl -X POST 'https://zenodo.org/api/deposit/depositions/$DEPOSITION_ID/actions/publish' \\"
echo "  -H 'Authorization: Bearer $ACCESS_TOKEN'"
echo ""
echo "========================================="
