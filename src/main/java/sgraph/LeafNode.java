package sgraph;

import com.jogamp.opengl.GL3;
import org.joml.Matrix4f;
import org.joml.Vector2f;
import org.joml.Vector3f;
import org.joml.Vector4f;
import raytracer.HitRecord;
import raytracer.Ray;
import util.Material;

import java.util.Map;
import java.util.Stack;

/**
 * This node represents the leaf of a scene graph. It is the only type of node that has
 * actual geometry to render.
 * @author Amit Shesh
 */
public class LeafNode extends AbstractNode
{
    /**
     * The name of the object instance that this leaf contains. All object instances are stored
     * in the scene graph itself, so that an instance can be reused in several leaves
     */
    protected String objInstanceName;
    /**
     * The material associated with the object instance at this leaf
     */
    protected util.Material material;

    protected String textureName;

    private static final float epsilon = 0.001f;

    public LeafNode(String instanceOf,IScenegraph graph,String name)
    {
        super(graph,name);
        this.objInstanceName = instanceOf;
    }



    /*
     *Set the material of each vertex in this object
     */
    @Override
    public void setMaterial(util.Material mat)
    {
        material = new util.Material(mat);
    }

    /**
     * Set texture ID of the texture to be used for this leaf
     * @param name
     */
    @Override
    public void setTextureName(String name)
    {
        textureName = name;
    }

    /*
     * gets the material
     */
    public util.Material getMaterial()
    {
        return material;
    }

    @Override
    public INode clone()
    {
        LeafNode newclone = new LeafNode(this.objInstanceName,scenegraph,name);
        newclone.setMaterial(this.getMaterial());
        return newclone;
    }


    /**
     * Delegates to the scene graph for rendering. This has two advantages:
     * <ul>
     *     <li>It keeps the leaf light.</li>
     *     <li>It abstracts the actual drawing to the specific implementation of the scene graph renderer</li>
     * </ul>
     * @param context the generic renderer context {@link sgraph.IScenegraphRenderer}
     * @param modelView the stack of modelview matrices
     * @throws IllegalArgumentException
     */
    @Override
    public void draw(IScenegraphRenderer context,Stack<Matrix4f> modelView) throws IllegalArgumentException
    {
        if (objInstanceName.length()>0)
        {
            context.drawMesh(objInstanceName,material,textureName,modelView.peek());
        }
    }


    /**
     * identify the object associated with the leaf (Box or Sphere)
     * and process the ray accordingly.
     */
    @Override
    public HitRecord calcIntersect(Ray ray, Stack<Matrix4f> modelView) {
        HitRecord hitRecord = new HitRecord();

        if (this.objInstanceName.indexOf("box") > -1) {
            hitRecord = this.boxIntersectionAlt(ray, modelView);
        } else if(this.objInstanceName.indexOf("sphere") > -1){
            hitRecord = this.sphereIntersection(ray, modelView);
        }else if(this.objInstanceName.indexOf("cylinder") > -1){
            hitRecord = this.cylinderIntersection(ray, modelView);
        }else{
            hitRecord = this.coneIntersection(ray, modelView);
        }

        return hitRecord;
    }


    /**
     *
     * @param ray ray from raycast
     * @param modelView to view transform
     * @return the hit record on box object
     */
    private HitRecord boxIntersection(Ray ray, Stack<Matrix4f> modelView) {
        HitRecord hitRecord = new HitRecord();
        Ray transRay = ray.rayTransform(new Matrix4f(modelView.peek()).invert());
        Vector4f startPoint = transRay.getStartPoint();
        Vector4f direction = transRay.getDirection();
        Vector4f normal = new Vector4f(0, 0, 0, 0);
        float minT = Float.NEGATIVE_INFINITY;
        float maxT = Float.POSITIVE_INFINITY;
        float temp = 0.0f, in = 0.0f, out = 0.0f;

        if (direction.x != 0) {
            float minTx = ((-0.5f - startPoint.x) / direction.x);
            float maxTx = ((0.5f - startPoint.x) / direction.x);
            if (minTx > maxTx) {
                temp = maxTx;
                maxTx = minTx;
                minTx = temp;
            }
            minT = Math.max(minT, Math.min(minTx, maxTx));
            maxT = Math.min(maxT, Math.max(minTx, maxTx));
        }

        if (transRay.getDirection().y != 0) {
            float minTy = ((-0.5f - startPoint.y) / direction.y);
            float maxTy = ((0.5f - startPoint.y) / direction.y);
            if (minTy > maxTy) {
                temp = maxTy;
                maxTy = minTy;
                minTy = temp;
            }
            minT = Math.max(minT, Math.min(minTy, maxTy));
            maxT = Math.min(maxT, Math.max(minTy, maxTy));
        }

        if (transRay.getDirection().z != 0) {
            float minTz = ((-0.5f - startPoint.z) / direction.z);
            float maxTz = ((0.5f - startPoint.z) / direction.z);
            if (minTz > maxTz) {
                temp = maxTz;
                maxTz = minTz;
                minTz = temp;
            }
            minT = Math.max(minT, Math.min(minTz, maxTz));
            maxT = Math.min(maxT, Math.max(minTz, maxTz));
        }

        //System.out.println(normal);

        // set normal
        if (maxT > 0 && minT <= maxT) {
            if (minT > 0) {
                in = minT;
            } else {
                in = maxT;
            }
            Vector4f intersectP = transRay.intersectObj(in);
            if (Math.abs(intersectP.x - 0.5f) < epsilon)
                normal.x = 1;
            else if (Math.abs(intersectP.x + 0.5f) < epsilon)
                normal.x = -1;
            else
                normal.x = 0;

            if (Math.abs(intersectP.y - 0.5f) < epsilon)
                normal.y = 1;
            else if (Math.abs(intersectP.y + 0.5f) < epsilon)
                normal.y = -1;
            else
                normal.y = 0;

            if (Math.abs(intersectP.z - 0.5f) < epsilon)
                normal.z = 1;
            else if (Math.abs(intersectP.z + 0.5f) < epsilon)
                normal.z = -1;
            else
                normal.z = 0;

            //System.out.println(Math.max(Math.abs(normal.x),Math.max(Math.abs(normal.y),Math.abs(normal.z))));
            Vector4f resNormal = new Matrix4f(modelView.peek()).invert().transpose().transform(normal);
            //System.out.println(resNormal);
            // 3D point of intersection in view coordinates
            Vector4f intersectP_view = ray.intersectObj(in);
            Vector4f intersectP_obj = transRay.intersectObj(in);


            hitRecord = new HitRecord(in, intersectP_view, intersectP_obj, resNormal, this.getMaterial());

            //set texture image coord
            hitRecord.setTextureName(this.textureName);
            if(Math.abs(Math.abs(intersectP_obj.x) - 0.5) < epsilon){
                if(intersectP_obj.x>0){
                    hitRecord.setTextureX((float)(1-(intersectP_obj.z+0.5)));
                    hitRecord.setTextureY((float)(intersectP_obj.y+0.5));
                }else{
                    hitRecord.setTextureX((float)(intersectP_obj.z+0.5));
                    hitRecord.setTextureY((float)(intersectP_obj.y+0.5));
                }

            }else if(Math.abs(Math.abs(intersectP_obj.y) - 0.5) < epsilon){
                if(intersectP_obj.y>0){
                    hitRecord.setTextureX((float)(1-(intersectP_obj.x+0.5)));
                    hitRecord.setTextureY((float)(intersectP_obj.z+0.5));
                }else {
                    hitRecord.setTextureX((float) (intersectP_obj.x + 0.5));
                    hitRecord.setTextureY((float) (intersectP_obj.z + 0.5));
                }
            }else {
                if(intersectP_obj.z>0){
                    hitRecord.setTextureX((float) ((intersectP_obj.x + 0.5)));
                    hitRecord.setTextureY((float) (intersectP_obj.y + 0.5));
                }else{
                    hitRecord.setTextureX((float) (1-(intersectP_obj.x + 0.5)));
                    hitRecord.setTextureY((float) (intersectP_obj.y + 0.5));
                }

            }
            if(hitRecord.getTextureX()>1)
                hitRecord.setTextureX(1f);
            if(hitRecord.getTextureY()>1)
                hitRecord.setTextureY(1);
        }
       // System.out.println("Box: " + hitRecord.isHit());
        return hitRecord;
    }


    /**
     * alternative way to calculate box hit, provided by professor
     * @param ray incoming ray
     * @param modelView view transform
     * @return hit record with box
     */
    private HitRecord boxIntersectionAlt(Ray ray,Stack<Matrix4f> modelView){
        HitRecord hitRecord = new HitRecord();
        Ray transRay = ray.rayTransform(new Matrix4f(modelView.peek()).invert());
        Vector4f startPoint = transRay.getStartPoint();
        Vector4f direction = transRay.getDirection();
        Vector4f normal = new Vector4f(0, 0, 0, 0);


        float tmaxX, tmaxY, tmaxZ;
        float tminX, tminY, tminZ;

        if (Math.abs(direction.x) < 0.0001f) {
            if ((startPoint.x > 0.5f) || (startPoint.x < -0.5f))
                return hitRecord;
            else {
                tminX = Float.NEGATIVE_INFINITY;
                tmaxX = Float.POSITIVE_INFINITY;
            }
        } else {
            float t1 = (-0.5f - startPoint.x) / direction.x;
            float t2 = (0.5f - startPoint.x) / direction.x;
            tminX = Math.min(t1, t2);
            tmaxX = Math.max(t1, t2);
        }

        if (Math.abs(direction.y) < 0.0001f) {
            if ((startPoint.y > 0.5f) || (startPoint.y < -0.5f)) {
                return hitRecord;
            } else {
                tminY = Float.NEGATIVE_INFINITY;
                tmaxY = Float.POSITIVE_INFINITY;
            }
        } else {
            float t1 = (-0.5f - startPoint.y) / direction.y;
            float t2 = (0.5f - startPoint.y) / direction.y;
            tminY = Math.min(t1, t2);
            tmaxY = Math.max(t1, t2);
        }

        if (Math.abs(direction.z) < 0.0001f) {
            if ((startPoint.z > 0.5f) || (startPoint.z < -0.5f)) {
                return hitRecord;
            } else {
                tminZ = Float.NEGATIVE_INFINITY;
                tmaxZ = Float.POSITIVE_INFINITY;
            }
        } else {
            float t1 = (-0.5f - startPoint.z) / direction.z;
            float t2 = (0.5f - startPoint.z) / direction.z;
            tminZ = Math.min(t1, t2);
            tmaxZ = Math.max(t1, t2);
        }

        float maxT, minT;
        minT = Math.max(tminX, Math.max(tminY, tminZ));
        maxT = Math.min(tmaxX, Math.min(tmaxY, tmaxZ));
        if (maxT > 0 && minT <= maxT) {
            float in;
            if (minT > 0) {
                in = minT;
            } else {
                in = maxT;
            }
            Vector4f intersectP = transRay.intersectObj(in);
            if (Math.abs(intersectP.x - 0.5f) < epsilon)
                normal.x = 1;
            else if (Math.abs(intersectP.x + 0.5f) < epsilon)
                normal.x = -1;
            else
                normal.x = 0;

            if (Math.abs(intersectP.y - 0.5f) < epsilon)
                normal.y = 1;
            else if (Math.abs(intersectP.y + 0.5f) < epsilon)
                normal.y = -1;
            else
                normal.y = 0;

            if (Math.abs(intersectP.z - 0.5f) < epsilon)
                normal.z = 1;
            else if (Math.abs(intersectP.z + 0.5f) < epsilon)
                normal.z = -1;
            else
                normal.z = 0;

            //System.out.println(Math.max(Math.abs(normal.x),Math.max(Math.abs(normal.y),Math.abs(normal.z))));
            Vector4f resNormal = new Matrix4f(modelView.peek()).invert().transpose().transform(normal);
            //System.out.println(resNormal);
            // 3D point of intersection in view coordinates
            Vector4f intersectP_view = ray.intersectObj(in);
            Vector4f intersectP_obj = transRay.intersectObj(in);


            hitRecord = new HitRecord(in, intersectP_view, intersectP_obj, resNormal, this.getMaterial());

            //set texture image coord
            hitRecord.setTextureName(this.textureName);
            if(Math.abs(Math.abs(intersectP_obj.x) - 0.5) < epsilon){
                if(intersectP_obj.x>0){
                    hitRecord.setTextureX((float)(1-(intersectP_obj.z+0.5)));
                    hitRecord.setTextureY((float)(intersectP_obj.y+0.5));
                }else{
                    hitRecord.setTextureX((float)(intersectP_obj.z+0.5));
                    hitRecord.setTextureY((float)(intersectP_obj.y+0.5));
                }

            }else if(Math.abs(Math.abs(intersectP_obj.y) - 0.5) < epsilon){
                if(intersectP_obj.y>0){
                    hitRecord.setTextureX((float)(1-(intersectP_obj.x+0.5)));
                    hitRecord.setTextureY((float)(intersectP_obj.z+0.5));
                }else {
                    hitRecord.setTextureX((float) (intersectP_obj.x + 0.5));
                    hitRecord.setTextureY((float) (intersectP_obj.z + 0.5));
                }
            }else {
                if(intersectP_obj.z>0){
                    hitRecord.setTextureX((float) ((intersectP_obj.x + 0.5)));
                    hitRecord.setTextureY((float) (intersectP_obj.y + 0.5));
                }else{
                    hitRecord.setTextureX((float) (1-(intersectP_obj.x + 0.5)));
                    hitRecord.setTextureY((float) (intersectP_obj.y + 0.5));
                }

            }
            if(hitRecord.getTextureX()>1)
                hitRecord.setTextureX(1f);
            if(hitRecord.getTextureY()>1)
                hitRecord.setTextureY(1);
        }
        // System.out.println("Box: " + hitRecord.isHit());
        return hitRecord;
    }

    /**
     *calculate sphere hit
     * @param ray ray from raycast
     * @param modelView to view transform
     * @return the hit record on sphere object
     */
    private HitRecord sphereIntersection(Ray ray, Stack<Matrix4f> modelView) {
        HitRecord hitRecord = new HitRecord();
        Ray transRay = ray.rayTransform(new Matrix4f(modelView.peek()).invert());

        Vector4f startPoint = transRay.getStartPoint();
        Vector4f direction = transRay.getDirection();
        float t1 = 0.0f, t2 = 0.0f;

        // (0, 0, 0) is center
        float argA = direction.x * direction.x + direction.y * direction.y + direction.z * direction.z;
        float argB = (direction.x * startPoint.x + direction.y * startPoint.y + direction.z * startPoint.z) * 2;
        float argC = startPoint.x * startPoint.x + startPoint.y * startPoint.y + startPoint.z * startPoint.z - 1;
        float sqrtVal = argB * argB - 4 * argA * argC;

        if (sqrtVal >= 0) {
            t1 = (float)((-1) * argB + (Math.sqrt(sqrtVal))) / (2 * argA);
            t2 = (float)((-1) * argB - (Math.sqrt(sqrtVal))) / (2 * argA);
            float minT = Math.min(t1, t2);
            if (minT <= 0) {
                minT = Math.max(t1, t2);
            }

            if (minT > 0) {
                //set normal
                Vector4f intersectP = transRay.intersectObj(minT);
                Vector4f normal = new Vector4f(intersectP);
                Matrix4f transpose = new Matrix4f(modelView.peek());
                transpose.invert().transpose().transform(normal);

                Vector4f intersectP_view = ray.intersectObj(minT);
                Vector4f intersectP_obj = transRay.intersectObj(minT);
                hitRecord = new HitRecord(minT, intersectP_view, intersectP_obj, normal, this.getMaterial());
                //set texture coordinate
                hitRecord.setTextureName(this.textureName);
                float phi = (float) Math.asin(intersectP_obj.y);
                phi = phi - (float) Math.PI / 2;
                float theta = (float) Math.atan2((double) intersectP.z, (double) intersectP_obj.x);
                theta  = theta +(float) Math.PI;
                float x = (float)(0.5-(theta / (2f * Math.PI))) ;
                float y = (1 - (phi / ((float) Math.PI)));
                hitRecord.setTextureX(x);
                hitRecord.setTextureY(y);
            }
        }
        //System.out.println("Sphere: " + hitRecord.isHit());
        return hitRecord;
    }


    /**
     * calculate cylinder intersection
     * @param ray incoming ray
     * @param modelView view transform
     * @return hit record
     */
    private HitRecord cylinderIntersection(Ray ray, Stack<Matrix4f> modelView){

        HitRecord capHitRecord = hitCylinderCap(ray, modelView);
        HitRecord wallHitRecord = hitCylinderWall(ray, modelView);


        if(capHitRecord.isHit() && wallHitRecord.isHit()){
            if(capHitRecord.getT()<wallHitRecord.getT())
                return capHitRecord;
            else
                return wallHitRecord;
        }else if(capHitRecord.isHit() && !wallHitRecord.isHit())
            return capHitRecord;
        else if(!capHitRecord.isHit() && wallHitRecord.isHit())
            return wallHitRecord;

        return new HitRecord();
    }

    /**
     * calculate if intersect with cylinder cap
     * @param ray incoming ray
     * @param modelView view transform
     * @return hit record
     */

    private HitRecord hitCylinderCap(Ray ray, Stack<Matrix4f> modelView){
        HitRecord hitRecord = new HitRecord();
        Ray transRay = ray.rayTransform(new Matrix4f(modelView.peek()).invert());

        Vector4f startPoint = transRay.getStartPoint();
        Vector4f direction = transRay.getDirection();
        Vector4f intersectionP = new Vector4f();
        float height = 1;
        float t = 0;

        if (transRay.getDirection().y != 0) {
            float bottomTy = ((0f - startPoint.y) / direction.y);
            float capTy = ((1.0f - startPoint.y) / direction.y);
            t = Math.min(bottomTy,capTy);
            if(Math.max(bottomTy, capTy)>epsilon){
                intersectionP = transRay.intersectObj(t);
                if(intersectionP.x*intersectionP.x+intersectionP.z*intersectionP.z > 1)
                    return hitRecord;
                Vector4f normal = new Vector4f();
                if(Math.abs(intersectionP.y) < epsilon)
                    normal = new Vector4f(0,-1,0,0);
                else
                    normal = new Vector4f(0,1,0,0);
                Matrix4f transpose = new Matrix4f(modelView.peek());
                transpose.invert().transpose().transform(normal);
                Vector4f intersectP_view = ray.intersectObj(t);
                Vector4f intersectP_obj = transRay.intersectObj(t);
                hitRecord = new HitRecord(t, intersectP_view, intersectP_obj, normal, this.getMaterial());
                hitRecord.setTextureName(this.textureName);
                float x = (intersectionP.x + 1)/2;
                float y = (intersectionP.z + 1)/2;
                if(x<0)
                    x = 0;
                if(x>1)
                    x = 1;
                if(y<0)
                    y = 0;
                if(y>1)
                    y =1;
                hitRecord.setTextureX(x);
                hitRecord.setTextureY(y);
                return  hitRecord;
            }else{
                return hitRecord;
            }
        }else{
            return hitRecord;
        }
    }

    /**
     * calculate if intersect with cylinder wall
     * @param ray incoming ray
     * @param modelView view transform
     * @return hit record
     */
    private HitRecord hitCylinderWall(Ray ray, Stack<Matrix4f> modelView){
        HitRecord hitRecord = new HitRecord();
        Ray transRay = ray.rayTransform(new Matrix4f(modelView.peek()).invert());

        Vector4f startPoint = transRay.getStartPoint();
        Vector4f direction = transRay.getDirection();
        Vector4f intersectionP = new Vector4f();
        float radius = 1;
        float t = 0;
        float a = direction.x*direction.x+direction.z*direction.z;
        float b = direction.x*startPoint.x +direction.z*startPoint.z;
        float c = startPoint.x*startPoint.x+startPoint.z*startPoint.z-radius*radius;

        float delta = b*b - a*c;
        if (delta < epsilon)
            return hitRecord;
        else{
            // nearest intersection
            t = (float)(-b - Math.sqrt (delta))/a;
            // t<0 means the intersection is behind the ray origin
            if (t<=epsilon)
                return hitRecord;
            else{
                intersectionP = transRay.intersectObj(t);
                if(intersectionP.y>1 || intersectionP.y<0)
                    return hitRecord;


                Vector4f normal = new Vector4f(intersectionP.x,0,intersectionP.z,0).normalize();
                Matrix4f transpose = new Matrix4f(modelView.peek());
                transpose.invert().transpose().transform(normal);
                Vector4f intersectP_view = ray.intersectObj(t);
                Vector4f intersectP_obj = transRay.intersectObj(t);
                hitRecord = new HitRecord(t, intersectP_view, intersectP_obj, normal, this.getMaterial());
                hitRecord.setTextureName(this.textureName);

                float theta = (float) Math.atan2((double) intersectionP.z, (double) intersectP_obj.x);
                theta  = theta +(float) Math.PI;
                float x = (float)(0.5-(theta / (2f * Math.PI))) ;
                hitRecord.setTextureX(x);
                hitRecord.setTextureY(intersectionP.y);


            }

        }
        return hitRecord;
    }


    /**
     * calculate cone intersection
     * @param ray incoming ray
     * @param modelView view transform
     * @return hit record
     */

    private  HitRecord coneIntersection(Ray ray, Stack<Matrix4f> modelView){
        HitRecord capHitRecord = hitConeCap(ray, modelView);
        HitRecord wallHitRecord = hitConeWall(ray, modelView);


        if(capHitRecord.isHit() && wallHitRecord.isHit()){
            if(capHitRecord.getT()<wallHitRecord.getT())
                return capHitRecord;
            else
                return wallHitRecord;
        }else if(capHitRecord.isHit() && !wallHitRecord.isHit())
            return capHitRecord;
        else if(!capHitRecord.isHit() && wallHitRecord.isHit())
            return wallHitRecord;

        return new HitRecord();
    }


    /**
     * calculate if hit cone wall
     * @param ray incoming ray
     * @param modelView view transform
     * @return hit record
     */
    private HitRecord hitConeWall(Ray ray,Stack<Matrix4f> modelView ){
        HitRecord hitRecord = new HitRecord();
        Ray transRay = ray.rayTransform(new Matrix4f(modelView.peek()).invert());

        Vector4f startPoint = transRay.getStartPoint();
        Vector4f direction = transRay.getDirection();
        Vector4f intersectionP = new Vector4f();
        float height = 1;
        float baseR = 1;
        float t = 0;

        Vector3f baseEdge = new Vector3f(baseR, 0, -height);
        Vector3f vertex = new Vector3f(0,0,0);
        float d = baseEdge.distance(vertex);
        float cosAlphaSq = height/d;
        float sinAlphaSq = baseR/d;

        float a = cosAlphaSq*direction.x*direction.x+cosAlphaSq*direction.z*direction.z-sinAlphaSq*direction.y*direction.y;
        float b = cosAlphaSq*direction.x*startPoint.x +cosAlphaSq*direction.z*startPoint.z-sinAlphaSq*direction.y*startPoint.y;
        float c = cosAlphaSq*startPoint.x*startPoint.x+cosAlphaSq*startPoint.z*startPoint.z-sinAlphaSq*startPoint.y*startPoint.y;

        double delta = b*b - a*c;
        // delta < 0 means no intersections
        if (delta < epsilon)
            return hitRecord;

        // nearest intersection
        t = (float)(-b - Math.sqrt (delta))/a;

        // t<0 means the intersection is behind the ray origin
        // which we don't want
        if (t<epsilon)
            return hitRecord;
        double y = startPoint.y+t*direction.y;
        // check if intersection point in cone range

        if (y < -height-epsilon || y > epsilon)
            return hitRecord;

        intersectionP = transRay.intersectObj(t);
        Vector3f v = new Vector3f(intersectionP.x,0,intersectionP.z).normalize();
        float nx = v.x*height/baseR;
        float ny = baseR/height;
        float nz = v.z*height/baseR;

        Vector4f normal = new Vector4f(nx,ny,nz,0).normalize();
        Matrix4f transpose = new Matrix4f(modelView.peek());
        transpose.invert().transpose().transform(normal);
        Vector4f intersectP_view = ray.intersectObj(t);
        Vector4f intersectP_obj = transRay.intersectObj(t);
        hitRecord = new HitRecord(t, intersectP_view, intersectP_obj, normal, this.getMaterial());
        hitRecord.setTextureName(this.textureName);
        float theta = (float) Math.atan2((double) intersectionP.z, (double) intersectP_obj.x);
        theta  = theta +(float) Math.PI;
        float x = (float)(0.5-(theta / (2f * Math.PI))) ;
        hitRecord.setTextureX(x);
        hitRecord.setTextureY(intersectionP.y);
        return  hitRecord;
    }


    /**
     * calculate cone cap/base
     * @param ray incoming ray
     * @param modelView view transform
     * @return hit record
     */

    private HitRecord hitConeCap(Ray ray, Stack<Matrix4f> modelView){
        HitRecord hitRecord = new HitRecord();
        Ray transRay = ray.rayTransform(new Matrix4f(modelView.peek()).invert());

        Vector4f startPoint = transRay.getStartPoint();
        Vector4f direction = transRay.getDirection();
        Vector4f intersectionP = new Vector4f();
        float height = 1;
        float t = 0;

        if (transRay.getDirection().y != 0) {
            t = ((-1.0f - startPoint.y) / direction.y);
            if(t>epsilon){
                intersectionP = transRay.intersectObj(t);
                if(intersectionP.x*intersectionP.x+intersectionP.z*intersectionP.z > 1)
                    return hitRecord;
                Vector4f normal = new Vector4f(0,-1,0,0);
                Matrix4f transpose = new Matrix4f(modelView.peek());
                transpose.invert().transpose().transform(normal);
                Vector4f intersectP_view = ray.intersectObj(t);
                Vector4f intersectP_obj = transRay.intersectObj(t);
                hitRecord = new HitRecord(t, intersectP_view, intersectP_obj, normal, this.getMaterial());
                hitRecord.setTextureName(this.textureName);
                float x = (intersectionP.x + 1)/2;
                float y = (intersectionP.z + 1)/2;
                if(x<0)
                    x = 0;
                if(x>1)
                    x = 1;
                if(y<0)
                    y = 0;
                if(y>1)
                    y =1;
                hitRecord.setTextureX(x);
                hitRecord.setTextureY(y);
                return  hitRecord;
            }else{
                return hitRecord;
            }
        }else{
            return hitRecord;
        }
    }



}
